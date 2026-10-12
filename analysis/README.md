# LUAD analysis: separate scripts and reusable functions

Open `Pathway_Modeling_Thesis.Rproj` or set your working directory to the repository root. Run the numbered scripts in order. Each script saves its result, so later scripts can be run separately without repeating earlier stages.

| Script | Purpose | Reusable function called |
|---|---|---|
| `01_prepare_luad.R` | Download, match, and prepare LUAD expression and survival data | `prepare_luad()` |
| `02_reactome_pathways.R` | Download Reactome sets and select 10 pathways | `load_reactome()` |
| `03_split_and_preprocess.R` | Split patients and preprocess expression using training patients | `split_patients()` |
| `04_pca_scores.R` | Calculate training PCA scores and project test patients | `pathway_pca()` |
| `05_gsva_scores.R` | Calculate GSVA scores with background genes retained | `pathway_gsva()` |
| `06_cox_models.R` | Fit Cox models for both scoring methods and predict test survival | `fit_survival_model(model="cox")` |
| `07_aft_models.R` | Fit Weibull AFT models for both scoring methods and predict test survival | `fit_survival_model(model="aft")` |
| `08_uno_cindex.R` | Calculate held-out Uno's C-index | `uno_cindex()` |
| `09_brier_score.R` | Calculate held-out Brier scores across time | `brier_score()` |
| `10_integrated_brier_score.R` | Calculate held-out integrated Brier scores | `integrated_brier_score()` |
| `11_compare_results.R` | Combine results and save the comparison table and plot | `compare_stage()` |

The Reactome script also saves the complete Reactome collection as **`data/pathways/reactome.gmt`**, alongside `kegg_medicus.gmt`. The GMT contains one pathway per line: name, description, and gene symbols separated by real tabs. The analysis selects its 10-pathway subset from this collection.

`00_setup.R` loads the shared functions and provides small file/score helpers. Your pathway functions remain in `pathway_functions/`, your metrics remain in `metric_functions/`, and your model function remains in `Model Function.R`. The numbered analysis scripts call these functions rather than containing duplicate implementations.

## In R / RStudio

Sourcing defines functions without starting a download or analysis. Call the stage function after sourcing its script:

```r
source("analysis/01_prepare_luad.R")
prepare_luad_stage()

source("analysis/02_reactome_pathways.R")
reactome_stage()

source("analysis/03_split_and_preprocess.R")
split_stage()

source("analysis/04_pca_scores.R")
pca_stage()

source("analysis/05_gsva_scores.R")
gsva_stage()

source("analysis/06_cox_models.R")
cox_stage()

source("analysis/07_aft_models.R")
aft_stage()

source("analysis/08_uno_cindex.R")
uno_cindex_stage()

source("analysis/09_brier_score.R")
brier_score_stage()

source("analysis/10_integrated_brier_score.R")
integrated_brier_score_stage()

source("analysis/11_compare_results.R")
compare_stage()
```

Inspect intermediate outputs in R, for example:

```r
readRDS("data/processed/luad_reactome_pilot/stages/reactome.rds")$sets
readRDS("data/processed/luad_reactome_pilot/stages/PCA_scores.rds")$scores$train
```

To run everything from R instead:

```r
source("analysis/luad_reactome.R")
prepare_luad_stage()
reactome_stage()
run_stages()
```

## From the cloud terminal

```bash
source /workspace/setup/activate.sh
export CURL_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
# Run one stage (after its prerequisite stages):
Rscript --vanilla analysis/02_reactome_pathways.R
# Or run the whole sequence:
Rscript --vanilla analysis/luad_reactome.R
```

The CA settings use the cloud machine's trusted HTTPS certificate bundle without disabling verification. Downloads require `linkedomics.org` and `zenodo.org`. Required packages are readr, msigdbr, GSVA, survival, and survAUC.

## Reusing your functions

`pathway_pca(expression_data, pathway_genes)` returns `patient_id` and `PC1`. Its `pca_fit` attribute stores the training fit; pass it with `pca_fit=` to score new patients without refitting. Expression data contain a `gene` column and numeric patient columns.

`pathway_gsva(expression_data, pathway_genes)` accepts either one gene vector or a named list of pathways. It retains the full background expression matrix.

`fit_survival_model(survival_data, pathway_scores, model)` accepts one or more score columns plus `patient_id`. It returns a prediction function. Call it with just new scores for the original linear predictions, `type="risk"` for risk ranking, or `times=..., type="survival"` for survival probabilities. The fitted model is in its `model` attribute.

The three metric functions retain their original positional arguments. Optional `training_time` and `training_event` supply a training censoring reference for held-out evaluation. Brier event weights use the censoring curve's left limit.

## Analysis choices and outputs

The pilot keeps seed 20261011, one event-stratified 70/30 split, and 10 reproducibly sampled Reactome pathways with 20–200 measured genes. Each model uses all 10 pathway scores together. Selection does not use survival associations or test performance. PCA transforms and expression imputation use training patients only. GSVA scores training and test cohorts separately; it is cohort-dependent and evaluates a held-out batch.

Metrics are Uno's C-index and Brier at the training follow-up 60th percentile, plus IBS over 51 points from its 10th to 60th percentiles. Both Cox and Weibull AFT supply survival probabilities. Times retain source units.

Final results are under `data/processed/luad_reactome_pilot/`, with reusable intermediate RDS files in `stages/`. Results include `performance.csv`, `selected_pathways.csv`, `patient_split.csv`, Brier curves and plot, fitted models, predictions, Cox diagnostics, checksums, and software versions.

Run `Rscript --vanilla analysis/luad_reactome.R --self-test` for the separate synthetic validation. It checks all four configurations, probability bounds, finite metrics, PCA projection, a hand-calculated Brier example, and censoring left limits. Synthetic numbers are not LUAD results.

This remains a single-split exploratory pilot, without confidence intervals, clinical-covariate adjustment, or external validation.
