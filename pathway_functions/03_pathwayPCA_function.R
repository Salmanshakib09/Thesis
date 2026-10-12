# expression_data: data frame with gene and one column per patient.
# Supply a fitted PCA object to project new patients without refitting.
pathway_pca <- function(expression_data, pathway_genes, pca_fit = NULL) {
  if (!"gene" %in% names(expression_data) || anyDuplicated(expression_data$gene)) stop("Supply unique gene identifiers")
  if (is.null(pca_fit)) {
    genes <- intersect(pathway_genes, expression_data$gene)
    if (length(genes) < 2) stop("Fewer than two genes present in the pathway")
  } else genes <- rownames(pca_fit$rotation)
  if (!all(genes %in% expression_data$gene)) stop("Missing genes needed for PCA projection")
  x <- t(as.matrix(expression_data[match(genes, expression_data$gene), setdiff(names(expression_data), "gene"), drop=FALSE]))
  colnames(x) <- genes
  if (any(!is.finite(x))) stop("PCA expression must be finite")
  if (is.null(pca_fit)) {
    pca_fit <- prcomp(x, center=TRUE, scale.=TRUE, rank.=1)
    scores <- pca_fit$x[, 1]
  } else scores <- predict(pca_fit, x)[, 1]
  result <- data.frame(patient_id=rownames(x), PC1=as.numeric(scores))
  attr(result, "pca_fit") <- pca_fit
  result
}
