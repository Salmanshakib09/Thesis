source("analysis/00_setup.R")
pca_stage <- function(directory=stage_dir) {
  x <- read_stage("split",directory)
  train <- test <- fits <- list()
  for (name in names(x$sets)) {
    tr <- pathway_pca(expression_frame(x$train_expression),x$sets[[name]])
    fits[[name]] <- attr(tr,"pca_fit")
    te <- pathway_pca(expression_frame(x$test_expression),x$sets[[name]],pca_fit=fits[[name]])
    train[[name]] <- tr$PC1; test[[name]] <- te$PC1
  }
  scores <- standardize_scores(do.call(cbind,train),do.call(cbind,test))
  save_stage(list(scores=scores,pca_fits=fits),"PCA_scores",directory)
}
if (sys.nframe()==0) pca_stage()
