# Repeated, nested validation of clinical, direct ridge, and stacked pathway models.
# Run from the project root after running 01_download_LUAD.R once.

library(tidyverse)
library(survival)
library(GSVA)

# Defaults use all eligible Reactome pathways. For a quick check, set
# REACTOME_MAX_PATHWAYS=50 before running Rscript. Zero means no pathway cap.
max_pathways <- as.integer(Sys.getenv("REACTOME_MAX_PATHWAYS", "0"))
n_repeats <- as.integer(Sys.getenv("REACTOME_REPEATS", "3"))
n_tuning_folds <- as.integer(Sys.getenv("REACTOME_TUNING_FOLDS", "3"))
n_stacking_folds <- as.integer(Sys.getenv("REACTOME_STACKING_FOLDS", "3"))
base_seed <- 20261011L
stacking_horizon <- 730 # Shared two-year prediction scale, fixed before validation.
penalties <- c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000)
output_folder <- Sys.getenv(
  "REACTOME_OUTPUT", "data/processed/reactome_prediction_comparison_scaled"
)
stopifnot(max_pathways >= 0, n_repeats >= 1,
          n_tuning_folds >= 2, n_stacking_folds >= 2)


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


# 2. Read pathway definitions and define training-only preparation ----

source("metric_functions/01_uno_cindex.R")
source("metric_functions/03_integrated_brier_score.R")

gmt_rows <- read_lines("data/pathways/reactome.gmt") |>
  discard(~ str_trim(.x) == "") |>
  str_split("\t")
stopifnot(length(gmt_rows) > 0, all(lengths(gmt_rows) >= 3))
pathway_genes <- map(gmt_rows, ~ unique(.x[-c(1, 2)])) |>
  set_names(map_chr(gmt_rows, 1))
stopifnot(!anyDuplicated(names(pathway_genes)))
pathway_genes <- pathway_genes[sort(names(pathway_genes))]
input_signature <- tools::md5sum(c(
  "data/processed/tcga_luad_processed.rds", clinical_file,
  "data/pathways/reactome.gmt", "Reactome_prediction_comparison.R"
))

make_folds <- function(data, k, seed) {
  set.seed(seed)
  folds <- integer(nrow(data))
  # Balance stage and event status so rare stages are represented internally.
  groups <- split(seq_len(nrow(data)), interaction(data$OS_event, data$stage, drop = TRUE))
  for (rows in groups) {
    order <- sample.int(k)
    assignments <- rep(order, length.out = length(rows))
    folds[rows] <- assignments[sample.int(length(rows))]
  }
  stopifnot(all(tabulate(folds, nbins = k) > 0))
  folds
}

clinical_design <- function(data) {
  model.matrix(~ age + sex + stage, data = data)[, -1, drop = FALSE]
}

# Convert risk to survival using matching centered references for numerical safety.
cox_survival <- function(model, newdata, times) {
  lp <- as.numeric(predict(model, newdata = newdata, type = "lp", reference = "sample"))
  hazard <- basehaz(model, centered = TRUE)
  h <- vapply(times, function(t) {
    values <- hazard$hazard[hazard$time <= t]
    if (length(values)) tail(values, 1) else 0
  }, numeric(1))
  probability <- exp(-exp(outer(lp, log(h), "+")))
  stopifnot(all(is.finite(probability)), all(probability >= 0 & probability <= 1))
  list(risk = lp, survival = probability)
}

# Compute both enrichment methods for every eligible pathway in one batch.
# GSVA uses each expression cohort separately, so this evaluates batch prediction.
prepare_pathways <- function(train, validation, horizon, label) {
  stopifnot(!any(train$patient_id %in% validation$patient_id))
  a <- expression_matrix[, train$patient_id, drop = FALSE]
  b <- expression_matrix[, validation$patient_id, drop = FALSE]
  keep <- apply(a, 1, sd) > 0
  a <- a[keep, , drop = FALSE]
  b <- b[keep, , drop = FALSE]
  genes <- map(pathway_genes, intersect, rownames(a))
  genes <- genes[between(lengths(genes), 20, 200)]
  if (max_pathways > 0) genes <- head(genes, max_pathways)
  stopifnot(length(genes) >= 2)

  score_batch <- function(mat) {
    list(
      gs = gsva(gsvaParam(mat, genes, kcdf = "Gaussian"),
                BPPARAM = BiocParallel::SerialParam(), verbose = FALSE),
      ss = gsva(ssgseaParam(mat, genes, normalize = FALSE),
                BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    )
  }
  scores_a <- score_batch(a)
  scores_b <- score_batch(b)
  failures <- list()
  models <- list()
  score_fits <- list()
  train_features <- list()
  valid_features <- list()
  risk <- list()

  for (pathway in names(genes)) {
    result <- tryCatch({
      g <- genes[[pathway]]
      pc <- prcomp(t(a[g, , drop = FALSE]), center = TRUE, scale. = TRUE, rank. = 1)
      pc_b <- predict(pc, newdata = t(b[g, , drop = FALSE]))[, 1]
      x <- cbind(PC1 = pc$x[, 1], GSVA = scores_a$gs[pathway, ], ssGSEA = scores_a$ss[pathway, ])
      z <- cbind(PC1 = pc_b, GSVA = scores_b$gs[pathway, ], ssGSEA = scores_b$ss[pathway, ])
      center <- colMeans(x)
      spread <- apply(x, 2, sd)
      stopifnot(all(is.finite(spread)), all(spread > 0))
      x <- scale(x, center = center, scale = spread)
      z <- scale(z, center = center, scale = spread)
      train_scores <- bind_cols(train, as_tibble(x))
      valid_scores <- bind_cols(validation, as_tibble(z))
      colnames(x) <- colnames(z) <- paste(pathway, colnames(x), sep = "|")
      train_features[[pathway]] <- x
      valid_features[[pathway]] <- z
      score_fits[[pathway]] <- list(pca = pc, center = center, scale = spread, genes = g)
      # Exclude failed/nonconverged first-stage models and record the reason.
      model <- withCallingHandlers(
        coxph(Surv(OS_time, OS_event) ~ PC1 + GSVA + ssGSEA + age + sex + stage,
              data = train_scores, x = TRUE, model = TRUE, na.action = na.fail),
        warning = function(w) stop(conditionMessage(w))
      )
      stopifnot(all(is.finite(coef(model))))
      environment(model$terms) <- globalenv()
      environment(model$formula) <- globalenv()
      attr(model$model, "terms") <- model$terms
      prediction <- cox_survival(model, valid_scores, horizon)
      # Log cumulative hazard gives fold models a common prediction scale.
      # Use log H_ref + lp directly to avoid underflow of survival probabilities.
      bh <- basehaz(model, centered = TRUE)
      h <- tail(bh$hazard[bh$time <= horizon], 1)
      stopifnot(length(h) == 1, is.finite(h), h > 0)
      pathway_risk <- log(h) + prediction$risk
      stopifnot(all(is.finite(pathway_risk)))
      list(model = model, pca = pc, center = center, scale = spread,
           genes = g, x = x, z = z, risk = pathway_risk)
    }, error = function(e) {
      failures[[pathway]] <<- tibble(context = label, pathway = pathway, reason = conditionMessage(e))
      NULL
    })
    if (!is.null(result)) {
      train_features[[pathway]] <- result$x
      valid_features[[pathway]] <- result$z
      risk[[pathway]] <- result$risk
      result$x <- result$z <- result$risk <- NULL
      models[[pathway]] <- result
    }
  }
  if (length(models) < 2) {
    print(list_rbind(failures) |> count(reason, sort = TRUE))
    stop("Fewer than two valid pathway models in ", label)
  }
  risks <- do.call(cbind, risk)
  colnames(risks) <- names(risk)
  message(label, ": ", length(models), " pathways; ", length(failures), " excluded")
  list(models = models, score_fits = score_fits, train_x = do.call(cbind, train_features),
       valid_x = do.call(cbind, valid_features), valid_risk = risks,
       direct_pathways = names(train_features), genes = genes[names(train_features)],
       failures = if (length(failures)) list_rbind(failures) else
         tibble(context = character(), pathway = character(), reason = character()))
}


# 3. Ridge Cox that supports more pathway scores than patients ----

# Standardize using training values, then rotate to the nonzero training span.
# An orthogonal rotation preserves the L2 penalty. Keeping this span gives the
# same minimum-norm ridge predictions without a p-by-p matrix inversion.
prepare_ridge <- function(x, z) {
  center <- colMeans(x)
  spread <- apply(x, 2, sd)
  keep <- is.finite(spread) & spread > 0
  stopifnot(sum(keep) >= 1)
  x <- scale(x[, keep, drop = FALSE], center[keep], spread[keep])
  z <- scale(z[, keep, drop = FALSE], center[keep], spread[keep])
  decomposition <- svd(x, nu = 0, nv = min(nrow(x), ncol(x)))
  rank <- sum(decomposition$d > max(decomposition$d) * 1e-8)
  stopifnot(rank > 0)
  rotation <- decomposition$v[, seq_len(rank), drop = FALSE]
  a <- x %*% rotation
  b <- z %*% rotation
  colnames(a) <- colnames(b) <- paste0("R", seq_len(rank))
  list(train = data.frame(scores = I(a)), validation = data.frame(scores = I(b)),
       columns = colnames(x), center = center[keep], scale = spread[keep],
       rotation = rotation)
}

fit_ridge <- function(design, train, penalty) {
  formula <- as.formula(paste0(
    "Surv(OS_time, OS_event) ~ ridge(scores, theta = ",
    penalty, ", scale = FALSE)"
  ), env = globalenv())
  model <- withCallingHandlers(
    coxph(formula, data = bind_cols(train, design$train),
          x = TRUE, model = TRUE, na.action = na.fail,
          control = coxph.control(iter.max = 100)),
    warning = function(w) stop(conditionMessage(w), call. = FALSE)
  )
  stopifnot(all(is.finite(coef(model))))
  model
}

# OOF scores are rebuilt from expression inside every stacking fold.
# This function is also called inside each tuning-training partition.
prepare_pipeline <- function(train, validation, horizon, seed, label) {
  # Reuse completed preparation only when data, code, settings, and patients match.
  cache_folder <- file.path(output_folder, "cache")
  dir.create(cache_folder, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_folder, paste0(gsub("/", "_", label), ".rds"))
  signature <- list(inputs = input_signature, train = train$patient_id,
                    validation = validation$patient_id, horizon = horizon,
                    seed = seed, max_pathways = max_pathways, folds = n_stacking_folds)
  if (file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (!is.null(cached) && identical(cached$signature, signature)) {
      message(label, ": reusing completed preparation")
      return(cached$pipeline)
    }
  }
  folds <- make_folds(train, n_stacking_folds, seed)
  oof <- list()
  failures <- list()
  for (fold in seq_len(n_stacking_folds)) {
    held <- which(folds == fold)
    run <- prepare_pathways(train[-held, ], train[held, ], horizon,
                            paste0(label, "/stack", fold))
    oof[[fold]] <- list(rows = held, risk = run$valid_risk)
    failures[[fold]] <- run$failures
  }
  full <- prepare_pathways(train, validation, horizon, paste0(label, "/refit"))
  common <- Reduce(intersect, c(map(oof, ~ colnames(.x$risk)), list(colnames(full$valid_risk))))
  stopifnot(length(common) >= 2)
  oof_risk <- matrix(NA_real_, nrow(train), length(common), dimnames = list(NULL, common))
  for (item in oof) oof_risk[item$rows, ] <- item$risk[, common, drop = FALSE]
  stopifnot(all(is.finite(oof_risk)))
  # Direct ridge uses all eligible scores, plus clinical predictors.
  direct <- prepare_ridge(
    cbind(full$train_x, clinical_design(train)),
    cbind(full$valid_x, clinical_design(validation))
  )
  stacked <- prepare_ridge(oof_risk, full$valid_risk[, common, drop = FALSE])
  pipeline <- list(direct = direct, stacked = stacked, full = full,
                   oof_risk = oof_risk, folds = folds, common_pathways = common,
                   failures = list_rbind(c(failures, list(full$failures))))
  saveRDS(list(signature = signature, pipeline = pipeline), cache_file)
  pipeline
}


# 4. Nested penalty tuning: tuning-validation outcomes never fit its pipeline ----

tune_penalties <- function(train, horizon, seed, label) {
  folds <- make_folds(train, n_tuning_folds, seed)
  tuning <- list()
  failure_tables <- list()
  for (fold in seq_len(n_tuning_folds)) {
    held <- which(folds == fold)
    inner <- train[-held, ]
    validation <- train[held, ]
    times <- seq(quantile(inner$OS_time, 0.1), quantile(inner$OS_time, 0.6), length.out = 21)
    pipeline <- prepare_pipeline(inner, validation, horizon, seed + fold,
                                 paste0(label, "/tune", fold))
    failure_tables[[fold]] <- pipeline$failures
    for (method in c("direct", "stacked")) {
      design <- pipeline[[method]]
      for (penalty in penalties) {
        # An unconverged candidate cannot be selected, even if its error looks low.
        loss <- tryCatch({
          model <- fit_ridge(design, inner, penalty)
          prediction <- cox_survival(model, design$validation, times)
          value <- integrated_brier_score(validation$OS_time, validation$OS_event,
                                           prediction$survival, times)
          list(value = value, error = NA_character_)
        }, error = function(e) list(value = Inf, error = conditionMessage(e)))
        tuning[[length(tuning) + 1L]] <- tibble(
          fold = fold, method = method, penalty = penalty,
          ibs = loss$value, fit_error = loss$error,
          n_validation = nrow(validation),
          n_direct_pathways = length(pipeline$full$direct_pathways),
          n_stacked_pathways = length(pipeline$common_pathways)
        )
      }
    }
  }
  scores <- list_rbind(tuning)
  selected <- scores |>
    group_by(method, penalty) |>
    summarise(mean_ibs = weighted.mean(ibs, n_validation),
              converged_folds = sum(is.finite(ibs)), .groups = "drop") |>
    filter(converged_folds == n_tuning_folds) |>
    arrange(method, mean_ibs, desc(penalty)) |>
    group_by(method) |>
    slice_head(n = 1) |>
    ungroup()
  stopifnot(setequal(selected$method, c("direct", "stacked")))
  list(selected = selected, scores = scores, folds = folds,
       failures = list_rbind(failure_tables))
}


# 5. Repeat the entire comparison on shared outer splits ----

evaluate <- function(model, newdata, train, test, times, label) {
  prediction <- cox_survival(model, newdata, times)
  probability <- prediction$survival
  stopifnot(all(apply(probability, 1, diff) <= 1e-10))
  curve <- map_dbl(seq_along(times), function(j) {
    brier_score(test$OS_time, test$OS_event, probability[, j], times[j])
  })
  list(
    metrics = tibble(
      model = label, n_training = nrow(train), n_test = nrow(test),
      test_events = sum(test$OS_event), time_start = first(times), time_end = last(times),
      uno_cindex = uno_cindex(test$OS_time, test$OS_event, prediction$risk,
                              last(times), train$OS_time, train$OS_event),
      brier_score = last(curve),
      integrated_brier_score = integrated_brier_score(test$OS_time, test$OS_event, probability, times)
    ),
    predictions = list(patient_id = test$patient_id, risk = prediction$risk, survival = probability),
    curve = tibble(model = label, time = times, brier_score = curve)
  )
}

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
excluded_patients |> write_csv(file.path(output_folder, "excluded_patients.csv"))
all_performance <- list()
all_tuning <- list()
all_selected <- list()
all_curves <- list()
all_differences <- list()

for (iteration in seq_len(n_repeats)) {
  # Shared event-stratified 70/30 split for all three methods.
  seed <- base_seed + iteration - 1L
  set.seed(seed)
  training_ids <- survival_data |>
    group_by(OS_event) |>
    slice_sample(prop = 0.7) |>
    ungroup() |>
    pull(patient_id)
  training_data <- filter(survival_data, patient_id %in% training_ids)
  test_data <- filter(survival_data, !patient_id %in% training_ids)
  times <- seq(quantile(training_data$OS_time, 0.1),
               quantile(training_data$OS_time, 0.6), length.out = 51)
  horizon <- stacking_horizon
  label <- paste0("repeat", iteration)
  tuning <- tune_penalties(training_data, horizon, seed + 100L, label)
  pipeline <- prepare_pipeline(training_data, test_data, horizon, seed + 200L, label)
  clinical_model <- coxph(Surv(OS_time, OS_event) ~ age + sex + stage,
                          data = training_data, x = TRUE, model = TRUE)
  results <- list(clinical = evaluate(clinical_model, test_data, training_data,
                                      test_data, times, "Clinical only"))
  ridge_models <- list()
  ridge_coefficients <- list()
  for (method in c("direct", "stacked")) {
    penalty <- tuning$selected$penalty[tuning$selected$method == method]
    model <- fit_ridge(pipeline[[method]], training_data, penalty)
    ridge_models[[method]] <- model
    design <- pipeline[[method]]
    weights <- as.numeric(design$rotation %*% coef(model))
    ridge_coefficients[[method]] <- tibble(
      method = method, predictor = design$columns,
      coefficient_per_training_sd = weights,
      coefficient_on_input_scale = weights / design$scale
    )
    results[[method]] <- evaluate(model, pipeline[[method]]$validation, training_data,
                                  test_data, times, paste("Ridge", method))
  }
  performance <- map(results, "metrics") |>
    list_rbind() |>
    mutate(iteration = iteration,
           n_direct_pathways = length(pipeline$full$direct_pathways),
           n_stacked_pathways = length(pipeline$common_pathways))
  stopifnot(all(is.finite(performance$uno_cindex)),
            all(is.finite(performance$integrated_brier_score)))
  differences <- performance |>
    filter(model != "Clinical only") |>
    transmute(iteration, model,
              delta_cindex = uno_cindex - results$clinical$metrics$uno_cindex,
              delta_ibs = integrated_brier_score - results$clinical$metrics$integrated_brier_score)
  curves <- map(results, "curve") |> list_rbind() |> mutate(iteration = iteration)
  folder <- file.path(output_folder, label)
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  performance |> write_csv(file.path(folder, "performance.csv"))
  differences |> write_csv(file.path(folder, "differences_from_clinical.csv"))
  list_rbind(ridge_coefficients) |> write_csv(file.path(folder, "ridge_coefficients.csv"))
  curves |> write_csv(file.path(folder, "brier_curves.csv"))
  tuning$scores |> write_csv(file.path(folder, "penalty_tuning.csv"))
  tuning$selected |> write_csv(file.path(folder, "selected_penalties.csv"))
  tuning$failures |> write_csv(file.path(folder, "excluded_tuning_pathway_models.csv"))
  pipeline$failures |> write_csv(file.path(folder, "excluded_pathway_models.csv"))
  tibble(pathway = pipeline$full$direct_pathways,
         measured_genes = lengths(pipeline$full$genes),
         first_stage_converged = pipeline$full$direct_pathways %in% names(pipeline$full$models),
         used_in_stack = pipeline$full$direct_pathways %in% pipeline$common_pathways) |>
    write_csv(file.path(folder, "selected_pathways.csv"))
  survival_data |>
    mutate(split = if_else(patient_id %in% training_ids, "training", "test")) |>
    write_csv(file.path(folder, "patient_split.csv"))
  bind_cols(training_data |> select(patient_id), tibble(fold = pipeline$folds),
            as_tibble(pipeline$oof_risk, .name_repair = "minimal")) |>
    write_csv(file.path(folder, "out_of_fold_predictions.csv"))
  bind_cols(test_data |> select(patient_id),
            as_tibble(pipeline$full$valid_risk, .name_repair = "minimal")) |>
    write_csv(file.path(folder, "test_pathway_predictions.csv"))
  # Keep transformations and fits needed to reproduce final predictions.
  list(clinical_model = clinical_model, ridge_models = ridge_models,
       ridge_coefficients = list_rbind(ridge_coefficients),
       pipeline = pipeline, results = results, tuning = tuning,
       training_data = training_data, test_data = test_data,
       times = times, seed = seed, settings = list(
         penalties = penalties, max_pathways = max_pathways,
         n_tuning_folds = n_tuning_folds, n_stacking_folds = n_stacking_folds,
         stacking_horizon = stacking_horizon, gsva_kernel = "Gaussian",
         ssgsea_normalize = FALSE, input_signature = input_signature
       )) |>
    saveRDS(file.path(folder, "analysis.rds"))
  all_performance[[iteration]] <- performance
  all_differences[[iteration]] <- differences
  all_curves[[iteration]] <- curves
  all_tuning[[iteration]] <- tuning$scores |> mutate(iteration = iteration)
  all_selected[[iteration]] <- tuning$selected |> mutate(iteration = iteration)
  message("Completed outer repeat ", iteration, "/", n_repeats)
  # Save progress after each outer split, so completed results survive interruption.
  list_rbind(all_performance) |> write_csv(file.path(output_folder, "performance.csv"))
}

performance <- list_rbind(all_performance)
performance_summary <- performance |>
  group_by(model) |>
  summarise(
    repeats = n(), mean_cindex = mean(uno_cindex), sd_cindex = sd(uno_cindex),
    mean_brier = mean(brier_score), sd_brier = sd(brier_score),
    mean_ibs = mean(integrated_brier_score), sd_ibs = sd(integrated_brier_score),
    .groups = "drop"
  )
performance_summary |> write_csv(file.path(output_folder, "performance_summary.csv"))
list_rbind(all_differences) |> write_csv(file.path(output_folder, "differences_from_clinical.csv"))
list_rbind(all_curves) |> write_csv(file.path(output_folder, "brier_curves.csv"))
list_rbind(all_tuning) |> write_csv(file.path(output_folder, "penalty_tuning.csv"))
list_rbind(all_selected) |> write_csv(file.path(output_folder, "selected_penalties.csv"))
capture.output(sessionInfo(), file = file.path(output_folder, "session_info.txt"))
performance_summary |> print(n = Inf, width = Inf)

# Repeated splits overlap: SD across repeats is descriptive, not a confidence interval.
# Each split uses its training 10th–60th follow-up percentiles for evaluation.
# Uno uses training censoring weights; Brier/IBS use evaluation-cohort weights.
# Full-cohort GSVA needs expression batches; this is not single-patient validation.
# View(performance); View(performance_summary); View(tuning$selected)
