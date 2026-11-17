# Five-fold validation of ten fixed Reactome pathway Cox models plus clinical Cox.
# Run from the project root after 01_download_LUAD.R. Does not fit stacked ridge.
# The pathways were fixed in the earlier exploratory pilot, not selected by CV outcomes.

library(tidyverse)
library(survival)
library(GSVA)
source("pathway_score_functions/03_pathwayPCA_function.R")
source("pathway_score_functions/04_GSVA_function.R")
source("metric_functions/01_uno_cindex.R")
source("metric_functions/02_brier_score.R")

k_folds <- 5L
fold_seed <- 20261117L
n_pathways <- 10L
# Prespecified pilot time window; identical prediction tasks in every fold.
evaluation_times <- seq(219.2, 809, length.out = 51)
output_folder <- "data/processed/reactome_covariates_cv"
figure_folder <- "figures/reactome_cv"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)

# 1. Assemble the same complete-case cohort as the original analysis ----

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


# 2. Fix the pathway identities before cross-validation ----
# Re-read membership from GMT; re-filter genes using each fold's training patients.
fixed_pathways <- c(
  "REACTOME_AQUAPORIN_MEDIATED_TRANSPORT",
  "REACTOME_CALNEXIN_CALRETICULIN_CYCLE",
  "REACTOME_CHROMATIN_ORGANIZATION",
  "REACTOME_CLASS_C_3_METABOTROPIC_GLUTAMATE_PHEROMONE_RECEPTORS",
  "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",
  "REACTOME_GOLGI_ASSOCIATED_VESICLE_BIOGENESIS",
  "REACTOME_NS1_MEDIATED_EFFECTS_ON_HOST_PATHWAYS",
  "REACTOME_RESPIRATORY_SYNCYTIAL_VIRUS_RSV_ATTACHMENT_AND_ENTRY",
  "REACTOME_SIGNALING_BY_FGFR_IN_DISEASE",
  "REACTOME_TELOMERE_C_STRAND_LAGGING_STRAND_SYNTHESIS"
)
gmt_rows <- read_lines("data/pathways/reactome.gmt") |>
  discard(~ str_trim(.x) == "") |>
  str_split("\t")
stopifnot(all(lengths(gmt_rows) >= 3))
all_pathways <- map(gmt_rows, ~ unique(.x[-c(1, 2)])) |>
  set_names(map_chr(gmt_rows, 1))
stopifnot(!anyDuplicated(names(all_pathways)),
          all(fixed_pathways %in% names(all_pathways)))
selected_genes <- all_pathways[fixed_pathways]

# Balance folds within each event/stage stratum without changing cohort order.
set.seed(fold_seed)
fold_ids <- integer(nrow(survival_data))
strata <- split(seq_len(nrow(survival_data)),
                interaction(survival_data$OS_event, survival_data$stage, drop = TRUE))
for (rows in strata) {
  assignments <- rep(sample.int(k_folds), length.out = length(rows))
  fold_ids[rows] <- assignments[sample.int(length(rows))]
}
stopifnot(all(tabulate(fold_ids, nbins = k_folds) > 0))
fold_assignments <- survival_data |> mutate(fold = fold_ids)
write_csv(fold_assignments, file.path(output_folder, "fold_assignments.csv"))
write_csv(excluded_patients, file.path(output_folder, "excluded_patients.csv"))

# 3. Fold-specific scoring and Cox fitting (same as the original first stage) ----

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

# 4. Common test metrics for each fold and model ----

predict_survival <- function(model, data) {
  risk <- as.numeric(predict(model, newdata = data, type = "lp", reference = "sample"))
  baseline <- basehaz(model, centered = TRUE)
  hazard <- map_dbl(evaluation_times, function(t) {
    tail(c(0, baseline$hazard[baseline$time <= t]), 1)
  })
  probabilities <- exp(-outer(exp(risk), hazard))
  stopifnot(all(is.finite(risk)), all(is.finite(probabilities)))
  list(risk = risk, survival = probabilities)
}

evaluate_model <- function(model, data, train, model_name, fold) {
  prediction <- predict_survival(model, data)
  censoring <- survfit(Surv(train$OS_time, 1 - train$OS_event) ~ 1)
  g_tau <- tail(c(1, censoring$surv[censoring$time <= max(evaluation_times)]), 1)
  if (g_tau <= 0 || sum(data$OS_time > max(evaluation_times)) == 0 ||
      sum(data$OS_event == 1 & data$OS_time <= max(evaluation_times)) < 2) {
    stop("Insufficient follow-up or censoring support in fold ", fold)
  }
  curve <- map_dbl(seq_along(evaluation_times), function(j) {
    brier_score(data$OS_time, data$OS_event,
                prediction$survival[, j], evaluation_times[j])
  })
  concordance <- uno_cindex(data$OS_time, data$OS_event, prediction$risk,
                            tau = max(evaluation_times),
                            training_time = train$OS_time,
                            training_event = train$OS_event)
  stopifnot(is.finite(concordance), all(is.finite(curve)))
  tibble(
    fold, model = model_name, train_n = nrow(train), test_n = nrow(data),
    test_deaths = sum(data$OS_event), uno_cindex = concordance,
    brier_score = tail(curve, 1),
    integrated_brier_score = sum(diff(evaluation_times) *
      (head(curve, -1) + tail(curve, -1)) / 2) / diff(range(evaluation_times))
  )
}

fold_results <- vector("list", k_folds)
for (fold in seq_len(k_folds)) {
  train <- survival_data[fold_ids != fold, ]
  test <- survival_data[fold_ids == fold, ]
  stopifnot(!any(train$patient_id %in% test$patient_id))
  message("Starting outer fold ", fold, "/", k_folds,
          ": train n = ", nrow(train), "; test n = ", nrow(test))
  models <- fit_pathways(train, test)
  clinical <- coxph(Surv(OS_time, OS_event) ~ age + sex + stage,
                    data = train, x = TRUE, model = TRUE, na.action = na.fail)
  stopifnot(all(is.finite(coef(clinical))))
  pathway_metrics <- imap(models, function(result, pathway) {
    stopifnot(identical(test$patient_id, result$test_scores$patient_id))
    data <- left_join(test, result$test_scores, by = "patient_id")
    evaluate_model(result$model, data, train, pathway, fold)
  }) |> list_rbind()
  fold_results[[fold]] <- bind_rows(
    evaluate_model(clinical, test, train, "Clinical-only Cox", fold), pathway_metrics
  )
  saveRDS(list(metrics = fold_results[[fold]], models = models,
               clinical_model = clinical, training_ids = train$patient_id,
               test_ids = test$patient_id),
           file.path(output_folder, paste0("fold_", fold, ".rds")))
  message("Completed outer fold ", fold, "/", k_folds)
}
metrics <- list_rbind(fold_results)
write_csv(metrics, file.path(output_folder, "fold_metrics.csv"))

# 5. Fold means and explicitly approximate standard errors ----
# SD/sqrt(k) ignores dependence from overlapping training sets. These bars are
# descriptive fold-based SEs, NOT confidence intervals or significance tests.
metric_labels <- c(uno_cindex = "Uno's C-index\n(higher is better)",
                   brier_score = "Brier score\n(lower is better)",
                   integrated_brier_score = "IBS\n(lower is better)")
long_metrics <- metrics |>
  pivot_longer(all_of(names(metric_labels)), names_to = "metric", values_to = "value")
summary_metrics <- long_metrics |>
  group_by(model, metric) |>
  summarise(mean = mean(value), sd = sd(value), n_folds = n(),
            se_naive = sd / sqrt(n_folds), .groups = "drop")
stopifnot(all(summary_metrics$n_folds == k_folds))
write_csv(summary_metrics, file.path(output_folder, "metric_summary.csv"))
write_csv(metrics, file.path(figure_folder, "fold_metrics.csv"))
write_csv(summary_metrics, file.path(figure_folder, "metric_summary.csv"))

# Pair the pathway and clinical results within the same fold for added value.
reference <- long_metrics |> filter(model == "Clinical-only Cox") |>
  select(fold, metric, reference = value)
differences <- long_metrics |> filter(model != "Clinical-only Cox") |>
  left_join(reference, by = c("fold", "metric")) |>
  mutate(improvement = if_else(metric == "uno_cindex", value - reference,
                                reference - value))
delta_summary <- differences |>
  group_by(model, metric) |>
  summarise(mean = mean(improvement), sd = sd(improvement), n_folds = n(),
            se_naive = sd / sqrt(n_folds), .groups = "drop")
write_csv(differences, file.path(output_folder, "paired_fold_improvements.csv"))
write_csv(delta_summary, file.path(output_folder, "paired_improvement_summary.csv"))
write_csv(delta_summary, file.path(figure_folder, "paired_improvement_summary.csv"))

pathway_label <- function(x) {
  x |> str_remove("^REACTOME_") |> str_replace_all("_", " ") |> str_wrap(30)
}
prepare_plot <- function(data) {
  data |> mutate(
    model = factor(model, levels = c(fixed_pathways, "Clinical-only Cox")),
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels)
  )
}
cv_plot <- prepare_plot(summary_metrics) |>
  ggplot(aes(mean, model, colour = model == "Clinical-only Cox")) +
  geom_segment(aes(x = mean - se_naive, xend = mean + se_naive, yend = model),
                linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_x", nrow = 1) +
  scale_y_discrete(labels = pathway_label) +
  scale_x_continuous(n.breaks = 3) +
  scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "firebrick"), guide = "none") +
  labs(title = "Five-fold held-out performance",
       subtitle = "Mean +/-1 naive fold SE",
       x = "Metric value", y = NULL) +
  theme_bw(base_size = 16) +
  theme(axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 12))

delta_plot <- prepare_plot(delta_summary) |>
  ggplot(aes(mean, model)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "firebrick") +
  geom_segment(aes(x = mean - se_naive, xend = mean + se_naive, yend = model),
                linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_x", nrow = 1) +
  scale_y_discrete(labels = pathway_label) +
  scale_x_continuous(n.breaks = 3) +
  labs(title = "Paired gain over clinical Cox",
       subtitle = "Mean +/-1 naive fold SE",
       x = "Improvement over clinical Cox", y = NULL) +
  theme_bw(base_size = 16) +
  theme(axis.text.y = element_text(size = 13),
        strip.text = element_text(size = 12))
for (extension in c("pdf", "png")) {
  ggsave(file.path(figure_folder, paste0("01_cv_performance.", extension)),
          cv_plot, width = 8.5, height = 9, dpi = 300)
  ggsave(file.path(figure_folder, paste0("02_cv_improvement.", extension)),
          delta_plot, width = 8.5, height = 9, dpi = 300)
}
list(folds = fold_assignments, metrics = metrics, summary = summary_metrics,
     paired_summary = delta_summary, times = evaluation_times,
     settings = list(k = k_folds, seed = fold_seed, pathways = fixed_pathways,
                     uncertainty = "Naive fold SD/sqrt(k); not a confidence interval")) |>
  saveRDS(file.path(output_folder, "cv_analysis.rds"))
capture.output(sessionInfo(), file = file.path(output_folder, "session_info.txt"))
message("Completed five-fold CV. Tables: ", output_folder, "; figures: ", figure_folder)
