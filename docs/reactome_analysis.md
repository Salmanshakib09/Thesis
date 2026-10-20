# Reactome analysis

Open the project in RStudio and run `01_download_LUAD.R` once if the processed
LUAD file is missing. Then source `Reactome_analysis.R` from the project root.
Required packages are tidyverse, survival, survAUC, and GSVA. GSVA is a
Bioconductor package. The GMT file is read directly; the pathwayPCA package is
not required because PCA uses base R's `prcomp()`.

## What the analysis fits

Each patient has ten pathway scores. All ten enter each of four joint models:
PCA–Cox, PCA–Weibull AFT, GSVA–Cox, and GSVA–Weibull AFT. There are four model
fits, rather than forty separate pathway fits. The event-stratified 70/30 split
is shared across the models. This remains a single-split exploratory analysis.

The settings at the top control the patient-split seed, pathway-selection seed,
number of pathways, GSVA kernel, and output folder. Both seeds initially retain
the previous value. Only these two steps draw random numbers in this workflow;
`PCA_simulation.R` has its own seed and is not sourced by this analysis.

## Preparation and scoring

Missing or infinite expression now stops the analysis with a preprocessing
message. No median imputation is performed. Genes with zero training variance
are removed, and pathways are eligible only when 20–200 retained genes remain.
Ten eligible pathways are sampled without using survival associations.

PCA is fitted once per pathway. Test patients use exactly the training gene
order, means, scales, and loadings. The existing `pathway_pca()` call still
returns a score table; `return_fit = TRUE` additionally returns the PCA fit.

GSVA computes all selected pathways together, estimating expression ranks only
once per cohort instead of repeating this work for each pathway. GSVA uses a
serial BiocParallel backend to keep resource use predictable. It retains
background genes and uses `kcdf = "Gaussian"` for both cohorts.
Previously, `auto` could choose direct ECDF for the larger training cohort and
Gaussian kernels for the smaller test cohort. This explicit choice can change
GSVA scores and model results. Training and test are still scored separately:
GSVA remains cohort-dependent and this is held-out batch evaluation, not a
validated scoring procedure for one future patient. The reusable GSVA helper
retains `auto` as its default for other existing callers.

Each pathway score is standardized using its training mean and standard
deviation. These values and the PCA fits are retained in `analysis.rds`.

## Survival probabilities and metrics

For Cox, the calculation is

`S(t | x) = exp(-H_ref(t) * exp(lp_centered))`.

`basehaz(fit, centered = TRUE)` supplies the cumulative hazard at the fitted
model's reference covariates. The default centered linear predictor from
`predict.coxph()` uses that same reference in these unstratified models.
Combining an uncentered hazard with a centered predictor would be incorrect.
The focused check compares these probabilities with `survfit(fit, newdata)`.
See the survival package's [basehaz documentation](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/basehaz.html)
and [predict.coxph documentation](https://stat.ethz.ch/R-manual/R-devel/library/survival/html/predict.coxph.html).

For Weibull AFT, with linear predictor `eta` and fitted scale `sigma`,
`S(t | x) = exp(-exp((log(t) - eta) / sigma))`. The check compares this with
`psurvreg()`. Risk for concordance is the Cox linear predictor or the negative
AFT linear predictor.

Uno's C-index now estimates censoring from training outcomes through survAUC's
`Surv.rsp` argument, evaluating concordance on test outcomes. Its old helper
signature still works, defaulting to the supplied evaluation outcomes when no
training outcomes are passed. See [UnoC documentation](https://cran.r-project.org/web/packages/survAUC/survAUC.pdf).
Brier and IBS retain test-cohort Kaplan–Meier censoring weights and left limits
at observed deaths, assuming independent censoring.

The time grid uses the 10th–60th percentiles of training follow-up. `start_time`
and `horizon` are included in the results: IBS averages over this interval,
not over zero to the horizon. Brier is reported at the horizon. RMSE is not
calculated. Small differences in one split do not establish a winning method.

## Inspect results

```r
performance |>
  View()

analysis_runs |>
  purrr::pluck("GSVA|cox", "model") |>
  summary()
```

Outputs remain under `data/processed/reactome_combined_analysis/` and are ignored
by Git. `analysis.rds` records models, scores, scaling values, retained genes,
times, and settings; `session_info.txt` records package versions.

Run the focused regression checks from the project root with:

```sh
Rscript tests/check_reactome_helpers.R
```
