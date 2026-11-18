# Five-fold validation of clinical, pathway-specific, and stacked ridge Cox models.
# Run from the project root after 01_download_LUAD.R. Then run the results script.
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
inner_folds <- 3L
stacking_folds <- 3L
ridge_penalties <- c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000)
stacking_horizon <- 730
# Two independent outer folds at a time on Mac/Linux; Windows uses one worker.
workers <- if (.Platform$OS.type == "windows") 1L else 2L
# Prespecified pilot time window; identical prediction tasks in every fold.
evaluation_times <- seq(219.2, 809, length.out = 51)
output_folder <- "data/processed/reactome_covariates_cv"
figure_folder <- "figures/reactome_cv"
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)
cache_folder <- file.path(output_folder, "nested_cache")
dir.create(cache_folder, recursive = TRUE, showWarnings = FALSE)

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


# 5. Leakage-safe nested stacking ----
# Changes to data, score functions or package versions invalidate cached fits.
cache_signature <- digest::digest(list(
  expression_matrix, survival_data, selected_genes,
  body(fit_pathways), body(fit_pathway_model), body(pathway_pca), body(pathway_gsva),
  packageVersion("GSVA"), packageVersion("survival"), packageVersion("stats")
))

cached_pathways <- function(train, validation, key) {
  file <- file.path(cache_folder, paste0(key, ".rds"))
  identity <- list(signature = cache_signature, training = train$patient_id,
                    validation = validation$patient_id)
  if (file.exists(file)) {
    saved <- readRDS(file)
    if (identical(saved$identity, identity)) return(saved$models)
  }
  models <- fit_pathways(train, validation)
  saveRDS(list(identity = identity, models = models), file)
  models
}

make_inner_folds <- function(data, n, seed) {
  set.seed(seed)
  ids <- integer(nrow(data))
  strata <- split(seq_len(nrow(data)), interaction(data$OS_event, data$stage, drop = TRUE))
  for (rows in strata) {
    values <- rep(sample.int(n), length.out = length(rows))
    ids[rows] <- values[sample.int(length(rows))]
  }
  stopifnot(all(tabulate(ids, nbins = n) > 0))
  ids
}

# Log cumulative hazard at a fixed horizon gives aligned risk features across
# fold models (centered linear predictors alone have arbitrary fold offsets).
stack_features <- function(models) {
  map(models, function(result) {
    baseline <- basehaz(result$model, centered = TRUE)
    h <- tail(baseline$hazard[baseline$time <= stacking_horizon], 1)
    stopifnot(length(h) == 1, is.finite(h), h > 0)
    log(h) + result$test_predictions$risk
  }) |> do.call(what = cbind)
}

prepare_stack <- function(train, validation, seed, key, models = NULL) {
  stopifnot(!any(train$patient_id %in% validation$patient_id))
  folds <- make_inner_folds(train, stacking_folds, seed)
  oof <- matrix(NA_real_, nrow(train), n_pathways,
                 dimnames = list(train$patient_id, fixed_pathways))
  for (fold in seq_len(stacking_folds)) {
    held <- which(folds == fold)
    base <- cached_pathways(train[-held, ], train[held, ],
                             paste0(key, "_oof", fold))
    oof[held, ] <- stack_features(base)
  }
  if (is.null(models)) {
    models <- cached_pathways(train, validation, paste0(key, "_validation"))
  }
  valid <- stack_features(models)
  center <- colMeans(oof)
  spread <- apply(oof, 2, sd)
  stopifnot(all(is.finite(oof)), all(is.finite(valid)), all(spread > 0))
  list(training = data.frame(scores = I(scale(oof, center, spread))),
       validation = data.frame(scores = I(scale(valid, center, spread))),
       oof = oof, folds = folds, center = center, scale = spread,
       training_ids = train$patient_id, validation_ids = validation$patient_id)
}

fit_ridge <- function(train, stack, penalty) {
  formula <- as.formula(paste0(
    "Surv(OS_time, OS_event) ~ ridge(scores, theta = ", penalty, ", scale = FALSE)"
  ), env = globalenv())
  model <- withCallingHandlers(
    coxph(formula, data = bind_cols(train, stack$training), model = TRUE, x = TRUE,
           na.action = na.fail, control = coxph.control(iter.max = 100)),
    warning = function(w) stop(conditionMessage(w), call. = FALSE)
  )
  stopifnot(all(is.finite(coef(model))))
  model
}

integrate_brier <- function(outcomes, prediction) {
  curve <- map_dbl(seq_along(evaluation_times), function(j) {
    brier_score(outcomes$OS_time, outcomes$OS_event,
                 prediction$survival[, j], evaluation_times[j])
  })
  sum(diff(evaluation_times) * (head(curve, -1) + tail(curve, -1)) / 2) /
    diff(range(evaluation_times))
}

run_outer_fold <- function(fold) {
  log_file <- file.path(output_folder, paste0("outer_", fold, ".log"))
  log_connection <- file(log_file, open = "wt")
  sink(log_connection, type = "message")
  on.exit({ sink(type = "message"); close(log_connection) })
  train <- survival_data[fold_ids != fold, ]
  test <- survival_data[fold_ids == fold, ]
  message("Outer fold ", fold, ": ", nrow(train), " training / ", nrow(test), " test")
  models <- cached_pathways(train, test, paste0("outer", fold))
  clinical <- coxph(Surv(OS_time, OS_event) ~ age + sex + stage,
                    data = train, x = TRUE, model = TRUE, na.action = na.fail)
  stopifnot(all(is.finite(coef(clinical))))

  # Each tuning-training set builds its OWN OOF features. Reusing an outer OOF
  # matrix here would allow tuning-validation outcomes into its predictors.
  tuning_folds <- make_inner_folds(train, inner_folds, fold_seed + 100L * fold)
  tuning <- map_dfr(seq_len(inner_folds), function(inner) {
    message("Outer ", fold, ": tuning fold ", inner, "/", inner_folds)
    held <- which(tuning_folds == inner)
    a <- train[-held, ]
    b <- train[held, ]
    stack <- prepare_stack(a, b, fold_seed + 1000L * fold + inner,
                             paste0("outer", fold, "_tune", inner))
    map_dfr(ridge_penalties, function(penalty) {
      loss <- tryCatch({
        model <- fit_ridge(a, stack, penalty)
        value <- integrate_brier(b, predict_survival(model, stack$validation))
        list(ibs = value, error = NA_character_)
      }, error = function(e) list(ibs = Inf, error = conditionMessage(e)))
      tibble(outer_fold = fold, inner_fold = inner, penalty,
              ibs = loss$ibs, n_validation = nrow(b), fit_error = loss$error)
    })
  })
  selected <- tuning |> group_by(penalty) |>
    summarise(mean_ibs = weighted.mean(ibs, n_validation), .groups = "drop") |>
    filter(is.finite(mean_ibs)) |> arrange(mean_ibs, desc(penalty)) |> slice_head(n = 1)
  stopifnot(nrow(selected) == 1)
  message("Outer ", fold, ": selected penalty ", selected$penalty)
  stack <- prepare_stack(train, test, fold_seed + 10000L * fold,
                           paste0("outer", fold, "_final"), models)
  ridge <- fit_ridge(train, stack, selected$penalty)

  pathway_metrics <- imap(models, function(result, pathway) {
    stopifnot(identical(test$patient_id, result$test_scores$patient_id))
    data <- left_join(test, result$test_scores, by = "patient_id")
    evaluate_model(result$model, data, train, pathway, fold)
  }) |> list_rbind()
  metrics <- bind_rows(
    evaluate_model(clinical, test, train, "Clinical-only Cox", fold),
    pathway_metrics,
    evaluate_model(ridge, bind_cols(test, stack$validation), train,
                     "Stacked ridge", fold)
  )
  # Retain held-out probability curves and patient predictions for audit.
  curve_rows <- map_dfr(c("Clinical-only Cox", "Stacked ridge"), function(name) {
    prediction <- if (name == "Clinical-only Cox") predict_survival(clinical, test) else
      predict_survival(ridge, stack$validation)
    tibble(fold, model = name, time = evaluation_times,
            brier_score = map_dbl(seq_along(evaluation_times), function(j) {
              brier_score(test$OS_time, test$OS_event,
                            prediction$survival[, j], evaluation_times[j])
            }))
  })
  predictions <- tibble(patient_id = test$patient_id, fold,
                         clinical_risk = predict_survival(clinical, test)$risk,
                         ridge_risk = predict_survival(ridge, stack$validation)$risk)
  saved <- list(metrics = metrics, models = models, clinical_model = clinical,
                  ridge_model = ridge, stack = stack, tuning = tuning,
                  tuning_folds = tuning_folds, selected_penalty = selected$penalty,
                  curves = curve_rows, predictions = predictions,
                  training_ids = train$patient_id, test_ids = test$patient_id)
  saveRDS(saved, file.path(output_folder, paste0("nested_fold_", fold, ".rds")))
  message("Completed outer fold ", fold)
  saved[c("metrics", "curves", "predictions", "tuning", "selected_penalty")]
}

# Forking only parallelizes independent OUTER folds. All GSVA calls remain serial.
fold_results <- parallel::mclapply(seq_len(k_folds), run_outer_fold,
                                  mc.cores = workers, mc.preschedule = FALSE)
if (any(vapply(fold_results, inherits, logical(1), "try-error"))) {
  stop("An outer fold failed. See outer_*.log; no incomplete summaries were generated.")
}
metrics <- map(fold_results, "metrics") |> list_rbind()
curves <- map(fold_results, "curves") |> list_rbind()
tuning <- map(fold_results, "tuning") |> list_rbind()
write_csv(metrics, file.path(output_folder, "fold_metrics.csv"))
write_csv(curves, file.path(output_folder, "fold_brier_curves.csv"))
write_csv(tuning, file.path(output_folder, "nested_ridge_tuning.csv"))
map(fold_results, "predictions") |> list_rbind() |>
  write_csv(file.path(output_folder, "held_out_predictions.csv"))
list(folds = fold_assignments, metrics = metrics, curves = curves, tuning = tuning,
       selected_penalties = map_dbl(fold_results, "selected_penalty"),
       times = evaluation_times,
       settings = list(k = k_folds, inner_folds = inner_folds,
                         stacking_folds = stacking_folds, seed = fold_seed,
                         pathways = fixed_pathways, penalties = ridge_penalties,
                         stacking_horizon = stacking_horizon,
                         cache_signature = cache_signature)) |>
  saveRDS(file.path(output_folder, "cv_analysis.rds"))
capture.output(sessionInfo(), file = file.path(output_folder, "session_info.txt"))
message("Completed nested five-fold CV. Run Reactome_covariates_cv_results.R next.")
