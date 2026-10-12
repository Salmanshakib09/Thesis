library(tidyverse)
library(GSVA)

pathway_gsva <- function(expression_data, pathway_genes) {
  genes_present <- pathway_genes |>
    intersect(expression_data$gene)
  if (length(genes_present) < 2) {
    stop("Fewer than two pathway genes are present in the expression data.")
  }
  
  expression_matrix <- expression_data |>
    filter(gene %in% genes_present) |>
    column_to_rownames("gene") |>
    as.matrix()
  
  pathway_list <- list(
    pathway = genes_present
  )
  
  gsva_parameters <- gsvaParam(
    expression_matrix,
    pathway_list
  )
  
  gsva_fit <- gsva_parameters |>
    gsva()
  
  pathway_scores <- gsva_fit |>
    t() |>
    as.data.frame() |>
    rownames_to_column("patient_id") |>
    rename(GSVA = pathway)
  
  return(pathway_scores)
}