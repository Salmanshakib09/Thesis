# Ten pathway-specific Cox models, followed by a stacked ridge Cox model.
# Run from the project root after running 01_download_LUAD.R once.

library(tidyverse)
library(survival)
library(GSVA)

"pathway_score_functions/03_pathwayPCA_function.R" |>
  source()
"pathway_score_functions/04_GSVA_function.R" |>
  source()
"pathway_score_functions/05_ssGSEA_function.R" |>
  source()

source("metric_functions/01_uno_cindex.R")
source("metric_functions/03_integrated_brier_score.R")

split_seed <- 20261011L
pathway_seed <- 20261011L
ridge_seed <- 20261021L
n_pathways <- 10L
n_folds <- 3L
ridge_penalties <- c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000)
stacking_horizon <- 730 # Common two-year prediction scale across fold models.
output_folder <- "data/processed/reactome_vector_covariates"


# 1. Load expression, survival outcomes, and clinical covariates ----

luad_data <- "data/processed/tcga_luad_processed.rds" |>
  readRDS() |>
  filter(
    is.finite(OS_time),
    OS_time > 0,
    OS_event %in% c(0, 1)
  )

# Cache the same clinical source used by 01_download_LUAD.R.
clinical_file <- "data/raw/luad_clinical_covariates.tsv"
clinical_url <- paste0(
  "https://linkedomics.org/data_download/TCGA-LUAD/",
  "Human__TCGA_LUAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi"
)

if (!file.exists(clinical_file)) {
  dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
  download.file(clinical_url, clinical_file, mode = "wb")
}

clinical_raw <- clinical_file |>
  read_tsv(
    col_types = cols(.default = col_character()),
    na = c("", "NA")
  )

clinical_matrix <- clinical_raw |>
  select(-1) |>
  as.matrix() |>
  t()

colnames(clinical_matrix) <- clinical_raw[[1]]

clinical_covariates <- clinical_matrix |>
  as_tibble(rownames = "patient_id") |>
  transmute(
    patient_id,
    age = as.numeric(years_to_birth),
    sex = factor(str_to_lower(str_trim(gender)), levels = c("female", "male")),
    stage = factor(
      str_to_lower(str_trim(pathologic_stage)),
      levels = c("stagei", "stageii", "stageiii", "stageiv")
    )
  )

stopifnot(!anyDuplicated(clinical_covariates$patient_id))

luad_data <- luad_data |>
  left_join(clinical_covariates, by = "patient_id")

# Record excluded patients; do not impute missing clinical values.
excluded_patients <- luad_data |>
  filter(!is.finite(age) | age <= 0 | is.na(sex) | is.na(stage)) |>
  select(patient_id, age, sex, stage)

luad_data <- luad_data |>
  filter(is.finite(age), age > 0, !is.na(sex), !is.na(stage))

survival_data <- luad_data |>
  select(patient_id, OS_time, OS_event, age, sex, stage)

expression_matrix <- luad_data |>
  select(-OS_time, -OS_event, -age, -sex, -stage) |>
  column_to_rownames("patient_id") |>
  as.matrix() |>
  t()

stopifnot(
  !anyNA(survival_data$patient_id),
  !anyDuplicated(survival_data$patient_id),
  !anyDuplicated(rownames(expression_matrix)),
  is.numeric(expression_matrix),
  all(is.finite(expression_matrix)),
  identical(colnames(expression_matrix), survival_data$patient_id)
)


# 2. Use one event-stratified training/test split for all pathways ----

set.seed(split_seed)

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

# Determine gene eligibility using training patients only.
gene_sd <- training_expression |>
  apply(MARGIN = 1, FUN = sd)
variable_genes <- names(gene_sd)[gene_sd > 0]

training_expression_data <- training_expression[variable_genes, , drop = FALSE] |>
  as_tibble(rownames = "gene", .name_repair = "minimal")
test_expression_data <- test_expression[variable_genes, , drop = FALSE] |>
  as_tibble(rownames = "gene", .name_repair = "minimal")


# 3. Select ten eligible Reactome pathways without outcome screening ----

gmt_rows <- "data/pathways/reactome.gmt" |>
  read_lines() |>
  discard(~ str_trim(.x) == "") |>
  str_split("\t")

stopifnot(length(gmt_rows) > 0, all(lengths(gmt_rows) >= 3))

pathway_genes <- gmt_rows |>
  map(~ unique(.x[-c(1, 2)])) |>
  set_names(map_chr(gmt_rows, 1)) |>
  map(intersect, variable_genes)

stopifnot(!anyDuplicated(names(pathway_genes)))

eligible_pathways <- tibble(
  pathway = names(pathway_genes),
  measured_genes = lengths(pathway_genes)
) |>
  filter(between(measured_genes, 20, 200)) |>
  arrange(pathway)

stopifnot(nrow(eligible_pathways) >= n_pathways)
set.seed(pathway_seed)

selected_pathways <- eligible_pathways |>
  slice_sample(n = n_pathways) |>
  arrange(pathway)
selected_genes <- pathway_genes[selected_pathways$pathway]


# 4. Reuse the same pathway fitting for full training data and internal folds ----

fit_pathway_model <- function(genes,
                              pathway,
                              training_data,
                              test_data,
                              training_expression_data,
                              test_expression_data,
                              gsva_training,
                              gsva_test,
                              ssgsea_training,
                              ssgsea_test) {

  # Fit PCA only on training patients, then project the test patients.
  pca_result <- training_expression_data |>
    pathway_pca(genes, return_fit = TRUE)
  pca_fit <- pca_result$fit

  pca_test <- pca_fit |>
    predict(
      newdata = test_expression_data |>
        column_to_rownames("gene") |>
        as.matrix() |>
        t()
    ) |>
    as.data.frame() |>
    rownames_to_column("patient_id") |>
    select(patient_id, PC1)

  training_scores <- pca_result$scores |>
    left_join(
      gsva_training |> select(patient_id, GSVA = all_of(pathway)),
      by = "patient_id"
    ) |>
    left_join(ssgsea_training |> select(patient_id, ssGSEA = all_of(pathway)),
              by = "patient_id")

  test_scores <- pca_test |>
    left_join(
      gsva_test |> select(patient_id, GSVA = all_of(pathway)),
      by = "patient_id"
    ) |>
    left_join(ssgsea_test |> select(patient_id, ssGSEA = all_of(pathway)),
              by = "patient_id")

  # Scale all three pathway scores using training estimates only.
  score_names <- c("PC1", "GSVA", "ssGSEA")
  score_center <- training_scores |>
    select(all_of(score_names)) |>
    map_dbl(mean)
  score_scale <- training_scores |>
    select(all_of(score_names)) |>
    map_dbl(sd)

  stopifnot(all(is.finite(score_scale)), all(score_scale > 0))

  training_scores <- training_scores |>
    mutate(
      across(
        all_of(score_names),
        ~ (.x - score_center[cur_column()]) / score_scale[cur_column()]
      )
    )
  test_scores <- test_scores |>
    mutate(
      across(
        all_of(score_names),
        ~ (.x - score_center[cur_column()]) / score_scale[cur_column()]
      )
    )

  model_training <- training_data |>
    left_join(training_scores, by = "patient_id")
  model_test <- test_data |>
    left_join(test_scores, by = "patient_id")

  # Each pathway has its own coefficients and baseline hazard. No penalty.
  fitted_model <- withCallingHandlers(
    model_training |>
    coxph(
      Surv(OS_time, OS_event) ~ PC1 + GSVA + ssGSEA + age + sex + stage,
      data = _,
      na.action = na.fail,
      model = TRUE,
      x = TRUE
    ),
    warning = function(w) stop(conditionMessage(w), call. = FALSE)
  )

  stopifnot(all(is.finite(coef(fitted_model))))
  # Keep the expression matrices out of saved formula environments.
  environment(fitted_model$terms) <- globalenv()
  environment(fitted_model$formula) <- globalenv()
  attr(fitted_model$model, "terms") <- fitted_model$terms

  test_risk <- fitted_model |>
    predict(newdata = model_test, type = "lp", reference = "sample") |>
    as.numeric()

  stopifnot(length(test_risk) == nrow(test_data), all(is.finite(test_risk)))

  list(
    model = fitted_model,
    pca_fit = pca_fit,
    score_center = score_center,
    score_scale = score_scale,
    training_scores = training_scores,
    test_scores = test_scores,
    test_predictions = tibble(patient_id = test_data$patient_id, risk = test_risk)
  )
}

fit_pathways <- function(train, validation) {
  # Refit gene filtering, PCA, and score scaling within each training partition.
  a <- expression_matrix[, train$patient_id, drop = FALSE]
  b <- expression_matrix[, validation$patient_id, drop = FALSE]
  keep <- apply(a, 1, sd) > 0
  a <- a[keep, , drop = FALSE]
  b <- b[keep, , drop = FALSE]
  genes <- map(selected_genes, intersect, rownames(a))
  stopifnot(length(genes) == n_pathways, all(lengths(genes) >= 2))
  a_data <- as_tibble(a, rownames = "gene", .name_repair = "minimal")
  b_data <- as_tibble(b, rownames = "gene", .name_repair = "minimal")
  gs_a <- pathway_gsva(a_data, genes, kcdf = "Gaussian")
  gs_b <- pathway_gsva(b_data, genes, kcdf = "Gaussian")

  # Batch ssGSEA is equivalent to the single-pathway helper with normalize = FALSE.
  score_ssgsea <- function(mat) {
    gsva(ssgseaParam(mat, genes, normalize = FALSE),
         BPPARAM = BiocParallel::SerialParam(), verbose = FALSE) |>
      t() |>
      as.data.frame() |>
      rownames_to_column("patient_id")
  }
  ss_a <- score_ssgsea(a)
  ss_b <- score_ssgsea(b)
  imap(genes, function(genes, pathway) {
    fit_pathway_model(genes, pathway, train, validation,
                      a_data, b_data, gs_a, gs_b, ss_a, ss_b)
  })
}

pathway_models <- fit_pathways(training_data, test_data)


# 5. Save the fitted models, associations, and held-out risk scores ----

pathway_test_predictions <- pathway_models |>
  imap(function(result, pathway) {
    result$test_predictions |>
      rename_with(~ pathway, -patient_id)
  }) |>
  reduce(left_join, by = "patient_id")

# These are mutually adjusted associations from training patients, not causal effects.
coefficient_table <- pathway_models |>
  imap(function(result, pathway) {
    estimates <- summary(result$model)$coefficients
    intervals <- summary(result$model)$conf.int
    tibble(
      pathway = pathway,
      term = rownames(estimates),
      coefficient = estimates[, "coef"],
      hazard_ratio = intervals[, "exp(coef)"],
      lower_95 = intervals[, "lower .95"],
      upper_95 = intervals[, "upper .95"],
      p_value = estimates[, "Pr(>|z|)"]
    )
  }) |>
  list_rbind()

patient_split <- survival_data |>
  mutate(split = if_else(patient_id %in% training_ids, "training", "test"))

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
selected_pathways |> write_csv(file.path(output_folder, "selected_pathways.csv"))
excluded_patients |> write_csv(file.path(output_folder, "excluded_patients.csv"))
patient_split |> write_csv(file.path(output_folder, "patient_split.csv"))
coefficient_table |> write_csv(file.path(output_folder, "coefficients.csv"))
pathway_test_predictions |> write_csv(file.path(output_folder, "test_risk_predictions.csv"))

list(
  models = pathway_models,
  selected_pathways = selected_pathways,
  patient_split = patient_split,
  excluded_patients = excluded_patients,
  settings = list(split_seed = split_seed, pathway_seed = pathway_seed, gsva_kernel = "Gaussian")
) |>
  saveRDS(file.path(output_folder, "analysis.rds"))

capture.output(sessionInfo(), file = file.path(output_folder, "session_info.txt"))

message(
  "Completed first stage: ", length(pathway_models), " pathway models; training n = ",
  nrow(training_data), "; test n = ", nrow(test_data),
  "; excluded for missing/invalid covariates = ", nrow(excluded_patients)
)

# Inspect one pathway's fitted model:
# pathway_models[[1]]$model |> summary()
# coefficient_table |> View()
# Score hazard ratios are per training SD; age is per year.
# Sex reference: female. Stage reference: stage I.
# Test scores are risk scores, not probabilities or predicted survival times.


# 6. Build out-of-fold pathway predictions for the ridge model ----

make_folds <- function(data, seed) {
  set.seed(seed)
  folds <- integer(nrow(data))
  groups <- split(seq_len(nrow(data)), interaction(data$OS_event, data$stage, drop = TRUE))
  for (rows in groups) {
    assignments <- rep(sample.int(n_folds), length.out = length(rows))
    folds[rows] <- assignments[sample.int(length(rows))]
  }
  stopifnot(all(tabulate(folds, nbins = n_folds) > 0))
  folds
}

# Unlike centered linear predictors alone, log cumulative hazards share a
# defined reference across folds. These are risk features, not probabilities.
pathway_predictions <- function(models) {
  map(models, function(result) {
    baseline <- basehaz(result$model, centered = TRUE)
    h <- tail(baseline$hazard[baseline$time <= stacking_horizon], 1)
    stopifnot(length(h) == 1, is.finite(h), h > 0)
    log(h) + result$test_predictions$risk
  }) |>
    do.call(what = cbind)
}

prepare_stack <- function(train, validation, seed, models = NULL) {
  stopifnot(!any(train$patient_id %in% validation$patient_id))
  folds <- make_folds(train, seed)
  oof <- matrix(NA_real_, nrow(train), n_pathways,
                dimnames = list(train$patient_id, names(selected_genes)))
  for (fold in seq_len(n_folds)) {
    held <- which(folds == fold)
    models_fold <- fit_pathways(train[-held, ], train[held, ])
    oof[held, ] <- pathway_predictions(models_fold)
    message("OOF fold ", fold, "/", n_folds, " completed")
  }
  if (is.null(models)) models <- fit_pathways(train, validation)
  valid <- pathway_predictions(models)
  center <- colMeans(oof)
  spread <- apply(oof, 2, sd)
  stopifnot(all(is.finite(oof)), all(is.finite(valid)),
            all(is.finite(spread)), all(spread > 0))
  list(training = data.frame(scores = I(scale(oof, center, spread))),
       validation = data.frame(scores = I(scale(valid, center, spread))),
       oof = oof, folds = folds, center = center, scale = spread)
}

fit_ridge <- function(train, stack, penalty) {
  formula <- as.formula(paste0(
    "Surv(OS_time, OS_event) ~ ridge(scores, theta = ", penalty, ", scale = FALSE)"
  ), env = globalenv())
  model <- withCallingHandlers(
    coxph(formula, data = bind_cols(train, stack$training),
          model = TRUE, x = TRUE, na.action = na.fail,
          control = coxph.control(iter.max = 100)),
    warning = function(w) stop(conditionMessage(w), call. = FALSE)
  )
  stopifnot(all(is.finite(coef(model))))
  model
}

ridge_predictions <- function(model, data, times) {
  risk <- as.numeric(predict(model, newdata = data, type = "lp", reference = "sample"))
  baseline <- basehaz(model, centered = TRUE)
  hazard <- map_dbl(times, function(t) {
    baseline |> filter(time <= t) |> pull(hazard) |> last(default = 0)
  })
  probability <- exp(-exp(outer(risk, log(hazard), "+")))
  stopifnot(all(is.finite(probability)),
            all(probability >= 0 & probability <= 1),
            all(apply(probability, 1, diff) <= 1e-10))
  list(risk = risk, survival = probability)
}


# 7. Tune the ridge penalty entirely within the training patients ----

# Rebuild the whole stack inside each tuning-training partition. Reusing one
# OOF table for tuning would let tuning-validation outcomes influence predictors.
tuning_folds <- make_folds(training_data, ridge_seed)
ridge_tuning <- map_dfr(seq_len(n_folds), function(fold) {
  message("Ridge tuning fold ", fold, "/", n_folds)
  held <- which(tuning_folds == fold)
  train <- training_data[-held, ]
  validation <- training_data[held, ]
  stack <- prepare_stack(train, validation, ridge_seed + fold)
  times <- seq(quantile(train$OS_time, 0.1), quantile(train$OS_time, 0.6), length.out = 21)
  map_dfr(ridge_penalties, function(penalty) {
    loss <- tryCatch({
      model <- fit_ridge(train, stack, penalty)
      prediction <- ridge_predictions(model, stack$validation, times)
      list(ibs = integrated_brier_score(validation$OS_time, validation$OS_event,
                                        prediction$survival, times), error = NA_character_)
    }, error = function(e) list(ibs = Inf, error = conditionMessage(e)))
    tibble(fold = fold, penalty = penalty, ibs = loss$ibs,
           n_validation = nrow(validation), fit_error = loss$error)
  })
})

selected_penalty <- ridge_tuning |>
  group_by(penalty) |>
  summarise(mean_ibs = weighted.mean(ibs, n_validation), .groups = "drop") |>
  filter(is.finite(mean_ibs)) |>
  arrange(mean_ibs, desc(penalty)) |>
  slice_head(n = 1)
stopifnot(nrow(selected_penalty) == 1)


# 8. Fit the final ten-pathway ridge model and evaluate untouched test patients ----

ridge_stack <- prepare_stack(training_data, test_data, ridge_seed + 1L,
                              models = pathway_models)
ridge_model <- fit_ridge(training_data, ridge_stack, selected_penalty$penalty)
evaluation_times <- seq(quantile(training_data$OS_time, 0.1),
                         quantile(training_data$OS_time, 0.6), length.out = 51)
ridge_test_predictions <- ridge_predictions(ridge_model, ridge_stack$validation,
                                             evaluation_times)
ridge_brier_curve <- tibble(
  time = evaluation_times,
  brier_score = map_dbl(seq_along(evaluation_times), function(j) {
    brier_score(test_data$OS_time, test_data$OS_event,
                ridge_test_predictions$survival[, j], evaluation_times[j])
  })
)
ridge_performance <- tibble(
  model = "Ten-pathway stacked ridge", n_pathways = n_pathways,
  n_training = nrow(training_data), n_test = nrow(test_data),
  test_events = sum(test_data$OS_event),
  time_start = first(evaluation_times), time_end = last(evaluation_times),
  uno_cindex = uno_cindex(test_data$OS_time, test_data$OS_event,
                          ridge_test_predictions$risk, last(evaluation_times),
                          training_data$OS_time, training_data$OS_event),
  brier_score = last(ridge_brier_curve$brier_score),
  integrated_brier_score = integrated_brier_score(
    test_data$OS_time, test_data$OS_event, ridge_test_predictions$survival, evaluation_times
  )
)
stopifnot(all(is.finite(unlist(ridge_performance[c(
  "uno_cindex", "brier_score", "integrated_brier_score"
)]))))
ridge_coefficients <- tibble(
  pathway = names(selected_genes),
  coefficient_per_training_sd = as.numeric(coef(ridge_model)),
  coefficient_on_input_scale = as.numeric(coef(ridge_model)) / ridge_stack$scale
)
ridge_performance |> write_csv(file.path(output_folder, "ridge_performance.csv"))
ridge_coefficients |> write_csv(file.path(output_folder, "ridge_coefficients.csv"))
ridge_brier_curve |> write_csv(file.path(output_folder, "ridge_brier_curve.csv"))
ridge_tuning |> write_csv(file.path(output_folder, "ridge_penalty_tuning.csv"))
selected_penalty |> write_csv(file.path(output_folder, "ridge_selected_penalty.csv"))
bind_cols(training_data |> select(patient_id), tibble(fold = ridge_stack$folds),
          as_tibble(ridge_stack$oof, .name_repair = "minimal")) |>
  write_csv(file.path(output_folder, "ridge_out_of_fold_predictions.csv"))
tibble(patient_id = test_data$patient_id, risk = ridge_test_predictions$risk) |>
  write_csv(file.path(output_folder, "ridge_test_risk_predictions.csv"))
list(model = ridge_model, stack = ridge_stack, performance = ridge_performance,
     coefficients = ridge_coefficients, predictions = ridge_test_predictions,
     tuning = ridge_tuning, tuning_folds = tuning_folds, selected_penalty = selected_penalty,
     times = evaluation_times, training_data = training_data, test_data = test_data,
     settings = list(n_folds = n_folds, penalties = ridge_penalties,
                     stacking_horizon = stacking_horizon, split_seed = split_seed,
                     ridge_seed = ridge_seed)) |>
  saveRDS(file.path(output_folder, "ridge_analysis.rds"))

ridge_performance |> print(width = Inf)
# View(ridge_performance); View(ridge_coefficients)
# Clinical covariates enter the pathway models; ridge combines their predictions.
# GSVA remains cohort-dependent: this evaluates batch prediction on one test split.
