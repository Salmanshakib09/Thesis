source("analysis/00_setup.R")
integrated_brier_score_stage <- function(directory=stage_dir) {
  x <- read_stage("split",directory); result <- list()
  for (key in c("PCA_Cox","PCA_AFT","GSVA_Cox","GSVA_AFT")) {
    pred <- read_stage(key,directory)
    result[[key]] <- integrated_brier_score(x$test$OS_time,x$test$OS_event,pred$survival,x$times,x$training$OS_time,x$training$OS_event)
    if (any(!is.finite(result[[key]]))) stop("Non-finite integrated_brier_score: ",key)
  }
  save_stage(result,"integrated_brier_score",directory)
}
if (sys.nframe()==0) integrated_brier_score_stage()
