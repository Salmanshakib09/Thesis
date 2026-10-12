# Convenience runner; each numbered script can also be run independently.
scripts <- c("01_prepare_luad","02_reactome_pathways","03_split_and_preprocess",
 "04_pca_scores","05_gsva_scores","06_cox_models","07_aft_models",
 "08_uno_cindex","09_brier_score","10_integrated_brier_score","11_compare_results")
for (script in scripts) source(paste0("analysis/",script,".R"))
run_stages <- function(directory=stage_dir) {
  split_stage(directory); pca_stage(directory); gsva_stage(directory)
  cox_stage(directory); aft_stage(directory)
  uno_cindex_stage(directory); brier_score_stage(directory); integrated_brier_score_stage(directory)
  compare_stage(directory)
}
if (sys.nframe()==0) {
  args <- commandArgs(trailingOnly=TRUE)
  if ("--self-test" %in% args) {
    source("analysis/test_modular_workflow.R")
    test_modular_workflow()
  } else {
    prepare_luad_stage()
    if (!"--prepare" %in% args) { reactome_stage(); run_stages() }
  }
}
