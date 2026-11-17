# Five-fold Reactome validation

From the project root, after preparing the TCGA LUAD data, run:

```r
source("Reactome_covariates_cross_validation.R")
```

This standalone script fits the ten fixed pilot pathway Cox models and the
clinical-only Cox reference in five event/stage-stratified outer folds. It does
not execute the original single-split script or its nested stacked-ridge stage.
It needs `data/processed/tcga_luad_processed.rds`, the clinical TSV (downloaded
and cached if absent), and `data/pathways/reactome.gmt`.

Each outer fold refits training gene filtering, PCA, score standardization,
and all Cox coefficients. GSVA uses separate training/test batches, retaining
its cohort dependence. ssGSEA uses `normalize = FALSE`. Pathway identities are
fixed from the exploratory pilot; they are not selected on CV test outcomes.
The previously chosen prediction window is fixed at 219.2–809 days, rather
than changing the prediction task from fold to fold. These prior pilot choices
make this exploratory internal validation, not independent confirmation.

The same folds are used for every model. Each patient is tested exactly once.
`fold_assignments.csv`, `fold_metrics.csv`, `metric_summary.csv`, paired
improvement tables, fitted fold models, and session information are saved to
`data/processed/reactome_covariates_cv/`. Figure PDFs and PNGs are written to
`figures/reactome_cv/`. The Results section references the two PDF figures.

Dots are equally weighted fold means. Bars are +/- one **naive fold-based SE**
(SD across five folds divided by sqrt(5)). Overlapping training sets induce
dependence, so these are descriptive approximate SEs, not confidence intervals
or significance tests. The clinical reference has its own uncertainty. The
second figure pairs each pathway metric with the clinical metric from the same
fold; positive differences favour pathways for all metrics.

Uno's C-index uses the relevant training-fold censoring distribution. Brier
scores use evaluation-fold censoring weights, matching the existing helper.
Evaluation fails clearly when follow-up/censoring support is insufficient.
The script stops on failed pathway-model fits rather than silently dropping
folds. Original single-split plots remain explicitly labelled in the report.

For a report update, commit the two plot PDFs and compact CSV summaries in
`figures/reactome_cv/`, then render `thesis.qmd`. These summaries are exported
automatically alongside the figures. Large RDS inputs and fitted models stay
untracked.
