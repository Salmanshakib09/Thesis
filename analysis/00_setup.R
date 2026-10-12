# Run/source scripts from the repository root. Sourcing never starts an analysis.
seed <- 20261011L
stage_dir <- "data/processed/luad_reactome_pilot/stages"
source("pathway_functions/03_pathwayPCA_function.R")
source("pathway_functions/04_GSVA_function.R")
source("Model Function.R")
source("metric_functions/01_uno_cindex.R")
source("metric_functions/02_brier_score.R")
source("metric_functions/03_integrated_brier_score.R")
save_stage <- function(x, name, directory=stage_dir) {
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
  saveRDS(x, file.path(directory, paste0(name, ".rds"))); invisible(x)
}
read_stage <- function(name, directory=stage_dir) {
  path <- file.path(directory, paste0(name, ".rds"))
  if (!file.exists(path)) stop("Missing ", name, ": run the preceding numbered script first")
  readRDS(path)
}
expression_frame <- function(x) data.frame(gene=rownames(x), x, check.names=FALSE)
standardize_scores <- function(train, test) {
  center <- colMeans(train); scale <- apply(train,2,sd)
  if (any(!is.finite(scale)) || any(scale<=0)) stop("Invalid score variance")
  list(train=sweep(sweep(train,2,center),2,scale,"/"),
       test=sweep(sweep(test,2,center),2,scale,"/"), center=center, scale=scale)
}
score_frame <- function(x, ids) data.frame(patient_id=ids, x, check.names=FALSE)
