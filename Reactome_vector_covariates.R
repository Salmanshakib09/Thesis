# One unregularized Cox model per pathway: PC1 + GSVA + ssGSEA + covariates.
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

split_seed <- 20261011L
pathway_seed <- 20261011L
n_pathways <- 10L
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


# 4. Calculate GSVA for all pathways, separately within each cohort ----

# Keep the same Gaussian kernel in both cohorts. GSVA is cohort-dependent.
gsva_training <- training_expression_data |>
  pathway_gsva(selected_genes, kcdf = "Gaussian")
gsva_test <- test_expression_data |>
  pathway_gsva(selected_genes, kcdf = "Gaussian")


# 5. Fit one Cox model using the three scores for each pathway ----

fit_pathway_model <- function(genes,
                              pathway) {

  # Fit PCA only on training patients, then project the test patients.
  pca_result <- training_expression_data |>
    pathway_pca(genes, return_fit = TRUE)
  pca_fit <- pca_result$fit

  pca_test <- pca_fit |>
    predict(
      newdata = t(test_expression[rownames(pca_fit$rotation), , drop = FALSE])
    ) |>
    as.data.frame() |>
    rownames_to_column("patient_id") |>
    select(patient_id, PC1)

  # The existing ssGSEA helper disables normalization across samples.
  ssgsea_training <- training_expression_data |>
    pathway_ssgsea(genes)
  ssgsea_test <- test_expression_data |>
    pathway_ssgsea(genes)

  training_scores <- pca_result$scores |>
    left_join(
      gsva_training |> select(patient_id, GSVA = all_of(pathway)),
      by = "patient_id"
    ) |>
    left_join(ssgsea_training, by = "patient_id")

  test_scores <- pca_test |>
    left_join(
      gsva_test |> select(patient_id, GSVA = all_of(pathway)),
      by = "patient_id"
    ) |>
    left_join(ssgsea_test, by = "patient_id")

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
  fitted_model <- model_training |>
    coxph(
      Surv(OS_time, OS_event) ~ PC1 + GSVA + ssGSEA + age + sex + stage,
      data = _,
      na.action = na.fail,
      model = TRUE,
      x = TRUE
    )

  stopifnot(all(is.finite(coef(fitted_model))))

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

pathway_models <- selected_genes |>
  imap(function(genes, pathway) {
    result <- fit_pathway_model(genes, pathway)
    message("Model completed: ", pathway)
    result
  })


# 6. Save the fitted models, associations, and held-out risk scores ----

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
  "Completed ", length(pathway_models), " pathway models; training n = ",
  nrow(training_data), "; test n = ", nrow(test_data),
  "; excluded for missing/invalid covariates = ", nrow(excluded_patients)
)

# Inspect one pathway's fitted model:
# pathway_models[[1]]$model |> summary()
# coefficient_table |> View()
# Score hazard ratios are per training SD; age is per year.
# Sex reference: female. Stage reference: stage I.
# Test scores are risk scores, not probabilities or predicted survival times.
# This is first-stage modeling only. Do not train a second-stage model on test
# predictions; generate out-of-fold training predictions for that later step.
