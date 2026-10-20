library(tidyverse)
library(GSVA)

pathway_ssgsea <- function(expression_data,
                           pathway_genes) {

  genes_present <- pathway_genes |>
    intersect(expression_data$gene)

  if (length(genes_present) < 2) {
    stop("Fewer than two pathway genes are present in the expression data.")
  }

  # Keep all measured genes as the background for ranking.
  expression_matrix <- expression_data |>
    column_to_rownames("gene") |>
    as.matrix()

  if (length(genes_present) >= nrow(expression_matrix)) {
    stop("Provide background genes outside the pathway.")
  }

  pathway_list <- list(
    pathway = genes_present
  )

  # Disable normalization across samples so scores do not depend on the batch.
  ssgsea_parameters <- ssgseaParam(
    expression_matrix,
    pathway_list,
    normalize = FALSE
  )

  pathway_scores <- ssgsea_parameters |>
    gsva(BPPARAM = BiocParallel::SerialParam()) |>
    t() |>
    as.data.frame() |>
    rownames_to_column("patient_id") |>
    rename(ssGSEA = pathway)

  return(pathway_scores)
}
