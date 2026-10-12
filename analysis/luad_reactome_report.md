# TCGA-LUAD: Reactome pathway prediction pilot

One complete pilot compared PCA and GSVA pathway scores with Cox and Weibull accelerated failure time (AFT) models. Each model used the same 10 Reactome pathways and the same held-out patients. PCA with Cox had the best observed values for all three metrics; differences were small and do not establish superiority.

## Cohort and reproducibility

LinkedOmics TCGA-LUAD clinical and log2 RSEM expression files supplied the data. Of 522 clinical records, 23 were excluded for missing/invalid survival endpoints or nonpositive follow-up, and 7 usable clinical records lacked expression data. The final cohort contained 492 patients, 178 deaths, and 19,988 measured genes.

An event-stratified 70/30 split with seed **20261011** assigned 343 patients to training (124 deaths) and 149 to testing (54 deaths). Patient identities were matched explicitly across clinical and expression columns. Exact membership is saved in the analysis outputs.

Reactome annotations came from **MSigDB 2026.1.Hs**, downloaded through `msigdbr`. Its archive matched the package's expected checksum. Ten pathways were sampled reproducibly from alphabetically ordered pathways containing 20–200 measured genes, using seed 20261011. Selection did not use survival associations or predictive performance. The subset is exploratory rather than a curated LUAD biological panel.

## Results on held-out patients

| Pathway scoring | Model | Uno's C-index ↑ | Brier score ↓ | Integrated Brier score ↓ |
|---|---|---:|---:|---:|
| PCA | Cox | **0.6455** | **0.2031** | **0.1387** |
| PCA | Weibull AFT | 0.6394 | 0.2037 | 0.1391 |
| GSVA | Cox | 0.6424 | 0.2086 | 0.1396 |
| GSVA | Weibull AFT | 0.6320 | 0.2099 | 0.1400 |

Uno's C-index and the single-time Brier score use a horizon of **807.4 in the source survival-time units**, chosen as the training follow-up 60th percentile. IBS averages Brier error over 51 equally spaced times from **187.4 to 807.4** (the training 10th to 60th percentiles). Calendar labels require confirmation of the source `overall_survival` time unit.

PCA–Cox exceeded GSVA–Cox by approximately 0.0031 in C-index and lowered IBS by approximately 0.0009. These descriptive differences are small. All C-indices exceeded 0.5, but a single test split without confidence intervals cannot establish statistical significance, generalizability, or a definitive winning configuration.

## Pathways

| Reactome pathway | Measured genes retained |
|---|---:|
| Aquaporin mediated transport | 52 |
| Ca dependent events | 36 |
| Chondroitin sulfate / dermatan sulfate metabolism | 37 |
| Class B/2 secretin family receptors | 90 |
| Co-inhibition by CTLA4 | 21 |
| Glutamate neurotransmitter release cycle | 23 |
| NR1H3/NR1H2 regulate gene expression linked to cholesterol transport and efflux | 33 |
| Respiratory syncytial virus attachment and entry | 23 |
| Signaling by FGFR2 | 71 |
| Telomere C-strand lagging-strand synthesis | 31 |

Exact database identifiers and feature mappings are saved in `selected_pathways.csv`. All 10 pathway scores enter each model jointly; this experiment fits four multivariable models rather than 40 single-pathway models.

## Methods and checks

- **PCA:** one first principal component per pathway. Gene centering, scaling, and PCA loadings were estimated only from training patients. Test patients were projected using those fitted transforms.
- **GSVA:** Gaussian-kernel GSVA with all nonconstant background genes retained. Training and test cohorts were scored separately so held-out expression did not influence training scores. This is cohort-dependent scoring and evaluates a held-out batch; it is not a validated single-patient scoring procedure.
- **Preprocessing:** training medians provide expression imputation; pathway-score centering and scaling use training patients only. No tuning or survival-based pathway selection was performed.
- **Models:** unpenalized Cox and Weibull AFT. Cox predictions use the training baseline cumulative hazard; AFT predictions use the fitted Weibull distribution. AFT linear predictions were negated when calculating risk-ranked C-index.
- **Metrics:** Uno's C-index uses the training censoring reference. Brier and IBS use inverse probability of censoring weights from training Kaplan–Meier censoring estimates, with left limits at event times. Both models supply survival probabilities for Brier metrics.
- **Validation:** four finite metric rows; disjoint patient sets; finite model coefficients; survival probabilities bounded by 0 and 1 and decreasing over time; centered training features; saved prediction dimensions. A separate synthetic four-model run and a hand-calculated Brier case passed.
- **Cox assumptions:** global Schoenfeld-residual tests gave p=0.2683 for PCA–Cox and p=0.7180 for GSVA–Cox. These tests did not detect a global departure from proportional hazards; they do not prove the assumption.

No model-convergence warnings were emitted. The first annotation download emitted a cross-filesystem rename warning; the package copied the verified archive successfully and completed normally.

## Deliverables and next research step

Run `analysis/luad_reactome.R` from the repository root; see `analysis/README.md` for environment activation and HTTPS trust configuration. Outputs are in `data/processed/luad_reactome_pilot/`: performance table, pathway list, patient split, Brier curves and PDF plot, fitted models and predictions, software versions, cohort accounting, input checksums, and Cox diagnostics.

The pilot fulfills the requested pathway × scoring × model × metric comparison. Its limits are one split, an exploratory pathway subset, cohort-dependent GSVA, no clinical-covariate adjustment, and no external validation or uncertainty estimates. The next research step is to discuss a biologically motivated pathway subset with the professor and use repeated or nested validation before drawing comparative conclusions.
