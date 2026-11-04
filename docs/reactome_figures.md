# Figures for the ten-pathway Reactome covariate analysis

Run the existing analysis first. To generate figures from its saved results,
run this in RStudio from the project root:

```r
"Reactome_covariates_figures.R" |>
  source()
```

This script does not change the existing analysis or rerun expression scoring.
It fits an additional clinical-only Cox reference using exactly the saved
training patients and evaluates it on the same training/test patients and time
grid. Figures and supporting CSV tables are saved in
`data/processed/reactome_vector_covariates/figures/`, in PNG and PDF formats.
The individual-pathway performance comparison opens in the Plots pane only when sourced
interactively. Use Files to open other exported figures.

| Figure | Research question | Context included | What it cannot establish |
| --- | --- | --- | --- |
| 01: Score associations | Which individual score coefficients are associated with mortality hazard conditional on the other predictors? | Training n/events; HR = 1 line; other scores and clinical adjustment; individual 95% CIs; BH adjustment across all 30 score tests; arbitrary PC1 sign. | Held-out prediction accuracy, causation, or joint added value of the three scores. |
| 02: Score correlations | Do PC1, GSVA and ssGSEA contain overlapping linear information? | Training-only Pearson correlations; -1 to +1 scale; r = 0 reference; arbitrary PC1 sign. | Association with survival, accuracy, or whether including a score improves prediction. |
| 03: Pathway versus clinical performance | Does each of the ten pathway-specific models improve prediction relative to clinical information alone on this split? | Clinical-only reference lines separately for each sample/metric; apparent training and held-out test rows; identical patients; event counts; Brier horizon and IBS interval; metric direction. | Statistical superiority, stable improvement across datasets, or the benefit of stacking. |
| 04: Final stack versus clinical | Does the combined ridge model improve on clinical-only prediction on this split? | Both models; training/test metrics; common cohort and time grid; clinical variables already enter first-stage models; training ridge uses OOF features but is still meta-model apparent performance. | Improvement over gene-level, direct pathway or random-pathway models, unless those comparisons are added. |
| 05: Brier curves | How does survival-probability error vary with follow-up? | Clinical reference curve; both training and test panels; follow-up units; censoring-weight convention; lower-is-better direction. | Calibration alone, a formal significance test, or performance beyond the plotted interval. |

Figures 04 and 05 require `ridge_analysis.rds` from the current analysis's final
ridge stage. If it is absent, the script clearly reports that they were not
generated. `figure_questions.csv` also records which figures were generated.
Never substitute a result from a different pathway collection or patient split.

## Interpretation rules

- Training metrics are apparent performance. Ridge's training features are OOF,
  but its training outcomes were still used to fit/tune the ridge model. These
  are not fully cross-validated metrics for the whole pipeline.
- Test metrics use untouched test outcomes for evaluation. Brier censoring
  weights are estimated from the evaluated sample, while Uno uses training
  censoring outcomes, following the existing metric helpers.
- The evaluation interval uses the training-follow-up 10th–60th percentiles,
  unless the saved ridge result supplies its grid. It does not start at zero.
- Each clinical reference includes age, sex and stage. No gene-expression data
  enter that baseline. All comparisons use the same complete-case cohort.
- One split provides point estimates, not confidence intervals. Do not select a
  pathway based on test performance and then present that same performance as
  independent validation.
- For the broader thesis, add gene-level models, biological-versus-random
  collections, scalar-versus-vector scores, and direct-versus-stacked models
  under shared validation. The current figures describe the ten-pathway pilot.
