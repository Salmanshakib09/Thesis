source("analysis/00_setup.R")
gsva_stage <- function(directory=stage_dir) {
  x <- read_stage("split",directory)
  # Separate cohort scoring retains all background genes in each cohort.
  train <- pathway_gsva(expression_frame(x$train_expression),x$sets)
  test <- pathway_gsva(expression_frame(x$test_expression),x$sets)
  scores <- standardize_scores(as.matrix(train[,names(x$sets),drop=FALSE]),as.matrix(test[,names(x$sets),drop=FALSE]))
  save_stage(list(scores=scores),"GSVA_scores",directory)
}
if (sys.nframe()==0) gsva_stage()
