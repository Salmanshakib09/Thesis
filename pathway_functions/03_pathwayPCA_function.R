library(tidyverse)
pathway_pca <- function(expression_data, pathway_genes) {
  genes_present<- pathway_genes |>
    intersect(expressio_data$gene)
  if(length(genes_present)<2) {
    stop("Fewer than two genes present in the pathway")
  }
  pca_matrix<- expression_data|>
    filter(gene %in% genes_present)|>
    column_to_rownames("gene")|>
    as.matrix()|>
    t()
  pca_fit<- pca_matrix|>
    prcomp(center=TRUE, scale.=TRUE)
  
  pathway_scores <- pca_fit$x |>
    as.data.frame() |>
    rownames_to_column("patient_id") |>
    select(patient_id, PC1)
  
  # Return pathway scores
  return(pathway_scores)
  
}