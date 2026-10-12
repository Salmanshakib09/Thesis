source("analysis/00_setup.R")
aft_stage <- function(directory=stage_dir) {
  x <- read_stage("split",directory)
  for (method in c("PCA","GSVA")) {
    scores <- read_stage(paste0(method,"_scores"),directory)$scores
    predict_model <- fit_survival_model(x$training,score_frame(scores$train,x$training$patient_id),model="aft")
    test <- score_frame(scores$test,x$test$patient_id)
    probability <- predict_model(test,times=x$times,type="survival")
    stopifnot(all(is.finite(probability)),all(probability>=0 & probability<=1),all(apply(probability,1,diff)<=1e-10))
    save_stage(list(model=attr(predict_model,"model"),risk=predict_model(test,type="risk"),
      survival=probability,patient_id=x$test$patient_id),paste(method,"AFT",sep="_"),directory)
  }
}
if (sys.nframe()==0) aft_stage()
