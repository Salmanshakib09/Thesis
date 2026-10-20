library(tidyverse)
library(GSVA)

pathway_gsva <- function(expression_data,
                         pathway_genes,
                         kcdf = "auto") {

  # A gene vector preserves the original patient_id + GSVA return format.
  # A named list scores several pathways against the same background at once.
  single_pathway <- !is.list(pathway_genes)
  pathway_list <- if (single_pathway) {
    list(GSVA = pathway_genes)
  } else {
    pathway_genes
  }

  if (is.null(names(pathway_list)) ||
      anyNA(names(pathway_list)) ||
      any(names(pathway_list) == "") ||
      anyDuplicated(names(pathway_list))) {
    stop("Provide unique names for the pathway list.")
  }

  pathway_list <- pathway_list |>
    map(intersect, expression_data$gene)

  if (length(pathway_list) == 0 || any(lengths(pathway_list) < 2)) {
    stop("Each pathway must contain at least two measured genes.")
  }

  expression_matrix <- expression_data |>
    column_to_rownames("gene") |>
    as.matrix()

  if (any(lengths(pathway_list) >= nrow(expression_matrix))) {
    stop("GSVA requires background genes outside each pathway.")
  }

  gsva_parameters <- gsvaParam(
    expression_matrix,
    pathway_list,
    kcdf = kcdf
  )

  pathway_scores <- gsva_parameters |>
    gsva(BPPARAM = BiocParallel::SerialParam()) |>
    t() |>
    as.data.frame() |>
    rownames_to_column("patient_id")

  return(pathway_scores)
}
