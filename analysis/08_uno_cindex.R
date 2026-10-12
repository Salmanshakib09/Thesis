source("analysis/00_setup.R")
uno_cindex_stage <- function(directory=stage_dir) {
  x <- read_stage("split",directory); result <- list()
  for (key in c("PCA_Cox","PCA_AFT","GSVA_Cox","GSVA_AFT")) {
    pred <- read_stage(key,directory)
    result[[key]] <- uno_cindex(x$test$OS_time,x$test$OS_event,pred$risk,tail(x$times,1),x$training$OS_time,x$training$OS_event)
    if (any(!is.finite(result[[key]]))) stop("Non-finite uno_cindex: ",key)
  }
  save_stage(result,"uno_cindex",directory)
}
if (sys.nframe()==0) uno_cindex_stage()
