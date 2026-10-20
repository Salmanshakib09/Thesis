# Run from the project root: Rscript tests/check_reactome_helpers.R
library(tidyverse)
library(survival)

source("pathway_score_functions/03_pathwayPCA_function.R")
source("Model Function.R")
source("metric_functions/01_uno_cindex.R")
source("metric_functions/03_integrated_brier_score.R")

# Shuffled pathway order must not create a second, differently oriented PCA.
set.seed(42)
expression_data <- matrix(rnorm(240), nrow = 6) |>
  as.data.frame() |>
  as_tibble() |>
  mutate(gene = paste0("gene", seq_len(6)), .before = 1)

genes <- c("gene5", "gene1", "gene3")
pca <- expression_data |>
  pathway_pca(genes, return_fit = TRUE)

legacy_scores <- expression_data |>
  pathway_pca(genes)

projection <- expression_data |>
  filter(gene %in% genes) |>
  column_to_rownames("gene") |>
  as.matrix() |>
  t() |>
  predict(pca$fit, newdata = _)

stopifnot(
  identical(legacy_scores, pca$scores),
  isTRUE(all.equal(unname(projection[, 1]), pca$scores$PC1))
)

# Batched GSVA must agree with the original single-pathway interface.
source("pathway_score_functions/04_GSVA_function.R")
sets <- list(first = genes, second = c("gene2", "gene4"))
batch <- pathway_gsva(expression_data, sets, kcdf = "Gaussian")
single <- pathway_gsva(expression_data, genes, kcdf = "Gaussian")
stopifnot(
  identical(batch$patient_id, single$patient_id),
  max(abs(batch$first - single$GSVA)) < 1e-10
)

# Check Cox survival probabilities independently using survfit(), including
# before the first event and beyond the last observed training time.
scores <- tibble(
  patient_id = paste0("patient", seq_len(80)),
  pathway_1 = rnorm(80),
  pathway_2 = rnorm(80)
)

outcomes <- scores |>
  transmute(
    patient_id,
    OS_time = pmax(1, round(rexp(80, exp(pathway_1 / 2)) * 100)),
    OS_event = rep(c(1, 1, 0, 1), 20)
  )

for (model in c("cox", "aft")) {
  predictor <- outcomes |>
    fit_survival_model(scores, model = model)

  fit <- attr(predictor, "model")
  new_scores <- scores |>
    slice(1:5)
  lp <- predictor(new_scores)
  times <- c(0.01, 20, 50, 100, max(outcomes$OS_time) + 10)

  if (model == "cox") {
    baseline <- basehaz(fit, centered = TRUE)
    hazard <- times |>
      map_dbl(function(t) {
        baseline |>
          filter(time <= t) |>
          pull(hazard) |>
          last(default = 0)
      })
    manual <- exp(-outer(exp(lp), hazard))
    reference <- survfit(fit, newdata = new_scores) |>
      summary(times = times, extend = TRUE)
    reference <- t(reference$surv)
  } else {
    manual <- exp(-exp(outer(-lp, log(times), "+") / fit$scale))
    reference <- times |>
      map(function(t) {
        1 - psurvreg(t, mean = lp, scale = fit$scale, distribution = "weibull")
      }) |>
      do.call(what = cbind, args = _)
  }

  stopifnot(max(abs(manual - reference)) < 1e-10)
}

# With a censoring time tied to a death, use G(t-) for that death.
# Contributions at t=2: .2^2 + .2^2 + 0 + .2^2/(2/3), divided by 4.
brier <- brier_score(
  actual_time = c(1, 2, 2, 3),
  actual_event = c(1, 1, 0, 1),
  predicted_survival = c(.2, .2, .8, .8),
  time_point = 2
)
stopifnot(abs(brier - .035) < 1e-12)

# Uno's wrapper must pass distinct training outcomes to survAUC.
index <- uno_cindex(
  outcomes$OS_time[61:80], outcomes$OS_event[61:80], scores$pathway_1[61:80],
  tau = 50,
  training_time = outcomes$OS_time[1:60],
  training_event = outcomes$OS_event[1:60]
)
expected <- survAUC::UnoC(
  Surv(outcomes$OS_time[1:60], outcomes$OS_event[1:60]),
  Surv(outcomes$OS_time[61:80], outcomes$OS_event[61:80]),
  scores$pathway_1[61:80], time = 50
)
stopifnot(isTRUE(all.equal(index, as.numeric(expected))))

# With no censoring and a constant prediction curve, IBS is that Brier score.
ibs <- integrated_brier_score(
  actual_time = c(1, 4),
  actual_event = c(1, 1),
  predicted_survival = matrix(.5, nrow = 2, ncol = 3),
  times = c(1, 2, 3)
)
stopifnot(abs(ibs - .25) < 1e-12)

bad_grid <- tryCatch(
  integrated_brier_score(c(1, 4), c(1, 1), matrix(.5, 2, 2), c(2, 1)),
  error = identity
)
stopifnot(inherits(bad_grid, "error"))

message("Passed: PCA projection, Cox/AFT probabilities, tied-time Brier, training-based Uno weights.")
