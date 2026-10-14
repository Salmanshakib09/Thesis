
# 1. Load packages and your functions ----

library(tidyverse)
library(survival)
library(GSVA)
library(pathwayPCA)

"pathway_score_functions/03_pathwayPCA_function.R" |>
  source()

"pathway_score_functions/04_GSVA_function.R" |>
  source()

source("Model Function.R")
source("metric_functions/01_uno_cindex.R")
source("metric_functions/02_brier_score.R")
source("metric_functions/03_integrated_brier_score.R")


# 2. Load and prepare LUAD data ----

# First run 01_download_LUAD.R to create this input file.
luad_data <- "data/processed/tcga_luad_processed.rds" |>
  readRDS() |>
  filter(
    is.finite(OS_time),
    OS_time > 0,
    OS_event %in% c(0, 1)
  )

survival_data <- luad_data |>
  select(patient_id, OS_time, OS_event)

# Expression matrix: genes in rows and patients in columns.
expression_matrix <- luad_data |>
  select(-OS_time, -OS_event) |>
  column_to_rownames("patient_id") |>
  as.matrix() |>
  t()

stopifnot(
  !anyDuplicated(survival_data$patient_id),
  !anyDuplicated(rownames(expression_matrix)),
  identical(colnames(expression_matrix), survival_data$patient_id)
)


# 3. Split patients into training and test groups ----

# Use the same event-balanced 70/30 split for all four combined models.
seed <- 20261011L
set.seed(seed)

training_ids <- survival_data |>
  group_by(OS_event) |>
  slice_sample(prop = 0.7) |>
  ungroup() |>
  pull(patient_id)

training_data <- survival_data |>
  filter(patient_id %in% training_ids)

test_data <- survival_data |>
  filter(!patient_id %in% training_ids)

training_expression <- expression_matrix[, training_data$patient_id, drop = FALSE]
test_expression <- expression_matrix[, test_data$patient_id, drop = FALSE]

# Learn missing-value replacements only from training patients.
gene_medians <- training_expression |>
  apply(
    MARGIN = 1,
    FUN = median,
    na.rm = TRUE
  )

usable_genes <- names(gene_medians)[is.finite(gene_medians)]
training_expression <- training_expression[usable_genes, , drop = FALSE]
test_expression <- test_expression[usable_genes, , drop = FALSE]

for (gene in usable_genes) {
  training_expression[gene, !is.finite(training_expression[gene, ])] <- gene_medians[gene]
  test_expression[gene, !is.finite(test_expression[gene, ])] <- gene_medians[gene]
}

# Remove genes with zero training variance before PCA or GSVA.
variable_genes <- training_expression |>
  apply(
    MARGIN = 1,
    FUN = sd
  )

variable_genes <- names(variable_genes)[variable_genes > 0]
training_expression <- training_expression[variable_genes, , drop = FALSE]
test_expression <- test_expression[variable_genes, , drop = FALSE]

# Your scoring functions expect a gene column and patient columns.
training_expression_data <- training_expression |>
  as_tibble(
    rownames = "gene",
    .name_repair = "minimal"
  )

test_expression_data <- test_expression |>
  as_tibble(
    rownames = "gene",
    .name_repair = "minimal"
  )


# 4. Read Reactome pathways and select a subset ----

reactome_pathways <- "data/pathways/reactome.gmt" |>
  read_gmt(
    description = TRUE
  )

pathway_genes <- reactome_pathways$pathways |>
  set_names(reactome_pathways$TERMS) |>
  map(
    intersect,
    rownames(expression_matrix)
  )

eligible_pathways <- tibble(
  pathway = names(pathway_genes),
  measured_genes = lengths(pathway_genes)
) |>
  filter(between(measured_genes, 20, 200)) |>
  arrange(pathway)

# Exploratory selection: it does not use survival associations or test results.
set.seed(seed)

selected_pathways <- eligible_pathways |>
  slice_sample(n = 10) |>
  arrange(pathway)

selected_genes <- pathway_genes[selected_pathways$pathway] |>
  map(
    intersect,
    variable_genes
  )

stopifnot(all(lengths(selected_genes) >= 2))


# 5. Calculate PCA and GSVA scores ----

# Standardize a score using its TRAINING mean and standard deviation.
standardize_scores <- function(training_scores,
                               test_scores,
                               score_name) {

  center <- training_scores[[score_name]] |>
    mean()

  spread <- training_scores[[score_name]] |>
    sd()

  stopifnot(is.finite(spread), spread > 0)

  list(
    training = training_scores |>
      mutate(across(all_of(score_name), ~ (.x - center) / spread)),
    test = test_scores |>
      mutate(across(all_of(score_name), ~ (.x - center) / spread))
  )
}

calculate_pathway_scores <- function(genes) {

  # Call your PCA helper for training scores.
  pca_training <- training_expression_data |>
    pathway_pca(
      pathway_genes = genes
    )

  # Your helper returns scores, not loadings. Retain a matching training PCA
  # fit separately so test patients can be projected without fitting test PCA.
  pca_fit <- training_expression[genes, , drop = FALSE] |>
    t() |>
    prcomp(
      center = TRUE,
      scale. = TRUE
    )

  stopifnot(isTRUE(all.equal(
    unname(pca_fit$x[, 1]),
    pca_training$PC1,
    tolerance = 1e-10
  )))

  pca_test <- pca_fit |>
    predict(
      newdata = t(test_expression[genes, , drop = FALSE])
    ) |>
    as.data.frame() |>
    rownames_to_column("patient_id") |>
    select(patient_id, PC1)

  # Call your GSVA helper separately for each cohort.
  # GSVA is cohort-dependent; this evaluates a held-out batch of patients.
  gsva_training <- training_expression_data |>
    pathway_gsva(
      pathway_genes = genes
    )

  gsva_test <- test_expression_data |>
    pathway_gsva(
      pathway_genes = genes
    )

  list(
    PCA = standardize_scores(pca_training, pca_test, "PC1"),
    GSVA = standardize_scores(gsva_training, gsva_test, "GSVA"),
    pca_fit = pca_fit
  )
}


# 6. Fit Cox/AFT and calculate all performance metrics ----

# Evaluation times are chosen from training follow-up, not test performance.
time_limits <- training_data$OS_time |>
  quantile(
    probs = c(0.1, 0.6)
  ) |>
  as.numeric()

evaluation_times <- seq(
  from = time_limits[1],
  to = time_limits[2],
  length.out = 51
)

evaluation_horizon <- evaluation_times |>
  last()

evaluate_model <- function(scoring,
                           model,
                           scores) {

  score_pair <- scores[[scoring]]

  stopifnot(
    identical(score_pair$training$patient_id, training_data$patient_id),
    identical(score_pair$test$patient_id, test_data$patient_id)
  )

  # Your model function returns a function that predicts linear values.
  prediction_function <- training_data |>
    fit_survival_model(
      pathway_scores = score_pair$training,
      model = model
    )

  linear_prediction <- score_pair$test |>
    prediction_function() |>
    as.numeric()

  # The returned function's environment contains its fitted survival model.
  fitted_model <- prediction_function |>
    environment() |>
    get(
      x = "fitted_model",
      envir = _
    )

  # Convert linear predictions to survival probabilities for Brier and IBS.
  if (model == "cox") {

    baseline_hazard <- fitted_model |>
      basehaz(
        centered = TRUE
      )

    baseline_at_times <- evaluation_times |>
      map_dbl(function(time_point) {
        baseline_hazard |>
          filter(time <= time_point) |>
          pull(hazard) |>
          last(default = 0)
      })

    survival_probability <- exp(
      -outer(exp(linear_prediction), baseline_at_times)
    )

    predicted_risk <- linear_prediction

  } else {

    # Weibull AFT: longer predicted survival corresponds to lower risk.
    survival_probability <- exp(
      -exp(
        outer(-linear_prediction, log(evaluation_times), "+") /
          fitted_model$scale
      )
    )

    predicted_risk <- -linear_prediction
  }

  stopifnot(
    all(is.finite(coef(fitted_model))),
    all(is.finite(survival_probability)),
    all(survival_probability >= 0 & survival_probability <= 1),
    all(apply(survival_probability, 1, diff) <= 1e-10)
  )

  # Call your metric functions as written. They estimate censoring on the
  # evaluated test patients. Survival times retain their source units.
  c_index <- test_data$OS_time |>
    uno_cindex(
      actual_event = test_data$OS_event,
      predicted_risk = predicted_risk,
      tau = evaluation_horizon
    )

  brier_scores <- evaluation_times |>
    seq_along() |>
    map_dbl(function(j) {
      test_data$OS_time |>
        brier_score(
          actual_event = test_data$OS_event,
          predicted_survival = survival_probability[, j],
          time_point = evaluation_times[j]
        )
    })

  integrated_brier <- test_data$OS_time |>
    integrated_brier_score(
      actual_event = test_data$OS_event,
      predicted_survival = survival_probability,
      times = evaluation_times
    )

  list(
    key = paste(scoring, model, sep = "|"),
    metrics = tibble(
      n_pathways = length(selected_genes),
      scoring = scoring,
      model = model,
      train_n = nrow(training_data),
      test_n = nrow(test_data),
      test_events = sum(test_data$OS_event),
      horizon = evaluation_horizon,
      uno_cindex = c_index,
      brier_score = last(brier_scores),
      integrated_brier_score = integrated_brier
    ),
    model = fitted_model,
    predictions = list(
      patient_id = test_data$patient_id,
      risk = predicted_risk,
      survival = survival_probability
    ),
    brier_curve = tibble(
      n_pathways = length(selected_genes),
      scoring = scoring,
      model = model,
      time = evaluation_times,
      brier_score = brier_scores
    )
  )
}


# 7. Combine the ten score columns and fit FOUR models ----

# Calculate scores for each pathway before fitting any survival model.
pathway_scores <- selected_genes |>
  imap(function(genes, pathway) {
    scores <- calculate_pathway_scores(genes)
    message("Scores completed: ", pathway)
    scores
  })

# Each resulting table has patient_id + TEN pathway-score columns.
combine_pathway_scores <- function(scoring,
                                   group) {

  pathway_scores |>
    imap(function(scores, pathway) {
      scores[[scoring]][[group]] |>
        rename_with(
          ~ pathway,
          -patient_id
        )
    }) |>
    reduce(
      left_join,
      by = "patient_id"
    )
}

combined_scores <- list(
  PCA = list(
    training = combine_pathway_scores("PCA", "training"),
    test = combine_pathway_scores("PCA", "test")
  ),
  GSVA = list(
    training = combine_pathway_scores("GSVA", "training"),
    test = combine_pathway_scores("GSVA", "test")
  )
)

stopifnot(
  ncol(combined_scores$PCA$training) == 11,
  ncol(combined_scores$GSVA$training) == 11
)

combinations <- crossing(
  scoring = c("PCA", "GSVA"),
  model = c("cox", "aft")
)

analysis_runs <- combinations |>
  pmap(function(scoring, model) {
    evaluate_model(scoring, model, combined_scores)
  })

analysis_runs <- analysis_runs |>
  set_names(map_chr(analysis_runs, "key"))

performance <- analysis_runs |>
  map("metrics") |>
  list_rbind() |>
  arrange(scoring, model)

stopifnot(
  nrow(performance) == 4,
  all(is.finite(performance$uno_cindex)),
  all(between(performance$uno_cindex, 0, 1)),
  all(is.finite(performance$brier_score)),
  all(is.finite(performance$integrated_brier_score))
)


# 8. Save the results and inspect them ----

output_folder <- "data/processed/reactome_combined_analysis"

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

performance |>
  write_csv(file.path(output_folder, "performance.csv"))

selected_pathways |>
  write_csv(file.path(output_folder, "selected_pathways.csv"))

survival_data |>
  mutate(
    split = if_else(patient_id %in% training_ids, "training", "test")
  ) |>
  write_csv(file.path(output_folder, "patient_split.csv"))

analysis_runs |>
  map("brier_curve") |>
  list_rbind() |>
  write_csv(file.path(output_folder, "brier_curves.csv"))

list(
  runs = analysis_runs,
  combined_scores = combined_scores,
  pathway_scores = pathway_scores,
  times = evaluation_times,
  seed = seed
) |>
  saveRDS(file.path(output_folder, "analysis.rds"))

capture.output(
  sessionInfo(),
  file = file.path(output_folder, "session_info.txt")
)

performance |>
  print(n = 4, width = Inf)

# In RStudio you can inspect the full results using:
# performance |> View()
# These are four joint ten-pathway models, not separate pathway models.
# This is a single-split exploratory study, not external validation.
