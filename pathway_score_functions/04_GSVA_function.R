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


# Score all pathways once per cohort, then scale with training estimates.
# GSVA distributions are still estimated separately in each cohort.
prepare_gsva_scores <- function(training_expression_data,
                                test_expression_data,
                                pathway_genes,
                                kcdf = "auto") {

  training_scores <- training_expression_data |>
    pathway_gsva(pathway_genes, kcdf = kcdf)

  test_scores <- test_expression_data |>
    pathway_gsva(pathway_genes, kcdf = kcdf)

  standardize_scores <- function(training_scores,
                                 test_scores,
                                 score_name) {

    center <- training_scores[[score_name]] |>
      mean()

    spread <- training_scores[[score_name]] |>
      sd()

    stopifnot(is.finite(spread), spread > 0)

    list(
      training = training_scores |>
        mutate(across(all_of(score_name), ~ (.x - center) / spread)),
      center = center,
      scale = spread,
      test = test_scores |>
        mutate(across(all_of(score_name), ~ (.x - center) / spread))
    )
  }

  names(training_scores)[-1] |>
    set_names() |>
    map(function(pathway) {
      training <- training_scores |>
        select(patient_id, GSVA = all_of(pathway))

      test <- test_scores |>
        select(patient_id, GSVA = all_of(pathway))

      standardize_scores(training, test, "GSVA")
    })
}
