# Reactome thesis analysis: five-fold nested validation

Run from the project root after preparing TCGA LUAD inputs:

```r
source("Reactome_covariates_cross_validation.R")
source("Reactome_covariates_cv_results.R")
```

The first script fits models; the second generates the three main figures,
cohort and performance tables, and an editable Results draft. Reporting does
not refit models. It does not overwrite your thesis prose when rerun.

## What is evaluated

The same five event/stage-stratified outer folds evaluate clinical-only Cox,
ten fixed Reactome pathway Cox models (PC1 + GSVA + ssGSEA + age + sex + stage),
and the final stacked ridge Cox model. Each patient is tested exactly once.
The pathways and 219.2–809-day window were fixed in the earlier pilot, making
this exploratory internal validation rather than independent confirmation.

For each outer training set:

1. Refit training gene filtering, PCA loadings, score scaling, and first-stage
   Cox coefficients. GSVA uses separate training/held-out batches; ssGSEA uses
   `normalize = FALSE`. GSVA remains cohort-dependent.
2. Select ridge strength using three inner tuning folds. In EACH tuning-training
   partition, build a separate three-fold OOF pathway-prediction matrix, fit
   ridge for the penalty candidates, and evaluate the tuning-validation patients.
   No single outer OOF matrix is reused for tuning.
3. Choose the penalty with minimum validation-size-weighted IBS; exact ties
   favour stronger regularization. Generate fresh three-fold OOF features for
   all outer training patients, fit ridge, and evaluate the untouched outer test
   fold. Features are log cumulative hazards at a common 730-day horizon, so
   their reference is aligned across first-stage fold models.

Clinical covariates enter the first-stage models; they are not added a second
time to the ridge stage. Inner fold assignments, candidate losses, selected
penalties, OOF features, models and test predictions are saved for audit.

The script stops if a first-stage fold fails or no candidate penalty can be
validated; it does not silently omit folds. First-stage cache files are keyed
by data, scoring code/package versions and patient partitions. Completed fits
can be reused when rerunning after interruption. Outer folds run two at a time
on Mac/Linux (`workers <- 1L` disables parallel work); Windows uses one worker.
Inner models are serial. Expect substantially longer runtime than first-stage
CV because every tuning partition must rebuild its stack.

## Outputs

- `data/processed/reactome_covariates_cv/`: raw model/checkpoint RDS files,
  fold assignments, held-out metrics/predictions, tuning records, and logs.
- `figures/reactome_cv/`: PDF/PNG main figures and compact metric summaries.
- `results/reactome_cv/`: patient table, performance table, pathway key,
  selected penalties, correlation/curve summaries, and `results_draft.qmd`.
- `thesis.qmd`: the reviewed Results narrative, tables and three figures.

The main figure sequence is performance, held-out Brier curves, then mean
absolute score correlations across outer training partitions. Paired gains
remain in the Results prose and CSV summaries, without a separate figure. Absolute correlations avoid cancellation when PC1 signs
reverse. These training correlations do not evaluate prediction accuracy.

The draft is a starting point for writing Results. Edit `thesis.qmd` for the
final narrative; rerunning the reporting script only updates the separate draft.
When incorporating new results, review the draft, copy the desired text into
the QMD, and render. Commit the figure PDFs and small CSV summaries with the
updated QMD/PDF. Large raw data and fitted model RDS files remain ignored.

## Uncertainty and interpretation

Dots are equally weighted outer-fold means. Bars/bands are +/- one naive
fold-based SE (fold SD/sqrt(5)). Overlapping training sets induce dependence,
so these are approximate descriptive SEs, NOT confidence intervals or formal
significance tests. Paired changes compare models on the same test fold.
Uno uses the outer training censoring distribution; Brier weights use the
corresponding evaluation fold, following the existing helpers. All models use
the same time grid. More methods, independent validation and more robust
uncertainty estimation remain future work; the current pilot does not answer
all gene-level, scalar/vector, and biological/random-pathway comparisons.

## Figure pagination

The thesis uses `fig-pos="H"` to keep the three main figures at their intended
positions, `needspace` to keep the Brier/correlation headings and discussion
with their plots, and a `placeins` float barrier before References. These LaTeX
packages are included in the QMD header. Generated Results drafts use the same
placement attributes; when copying a draft to another QMD, keep those header
packages and the barrier before its bibliography. The reporting script no
longer creates the paired-gain plot, but still exports its numerical summary.
