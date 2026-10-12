source("analysis/00_setup.R")
brier_score_stage <- function(directory=stage_dir) {
  x <- read_stage("split",directory); result <- list()
  for (key in c("PCA_Cox","PCA_AFT","GSVA_Cox","GSVA_AFT")) {
    pred <- read_stage(key,directory)
    result[[key]] <- vapply(seq_along(x$times),function(j) brier_score(x$test$OS_time,x$test$OS_event,pred$survival[,j],x$times[j],x$training$OS_time,x$training$OS_event),numeric(1))
    if (any(!is.finite(result[[key]]))) stop("Non-finite brier_score: ",key)
  }
  save_stage(result,"brier_score",directory)
}
if (sys.nframe()==0) brier_score_stage()
