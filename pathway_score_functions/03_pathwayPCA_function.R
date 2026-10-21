library(tidyverse)

pathway_pca <- function(expression_data,
                        pathway_genes,
                        return_fit = FALSE) {

  genes_present <- pathway_genes |>
    intersect(expression_data$gene)

  if (length(genes_present) < 2) {
    stop("Fewer than two genes present in the pathway.")
  }

  pca_matrix <- expression_data |>
    filter(gene %in% genes_present) |>
    column_to_rownames("gene") |>
    as.matrix() |>
    t()

  pca_fit <- pca_matrix |>
    prcomp(center = TRUE, scale. = TRUE)

  pathway_scores <- pca_fit$x |>
    as.data.frame() |>
    rownames_to_column("patient_id") |>
    select(patient_id, PC1)

  # Existing calls still return a score table.
  if (return_fit) {
    return(list(scores = pathway_scores, fit = pca_fit))
  }

  return(pathway_scores)
}


# Prepare training and test PC1 scores using one training PCA fit.
prepare_pca_scores <- function(training_expression_data,
                               test_expression_data,
                               pathway_genes) {

  pca_result <- training_expression_data |>
    pathway_pca(pathway_genes, return_fit = TRUE)

  pca_fit <- pca_result$fit
  pca_training <- pca_result$scores

  test_matrix <- test_expression_data |>
    column_to_rownames("gene") |>
    as.matrix()

  pca_test <- pca_fit |>
    predict(
      newdata = t(test_matrix[rownames(pca_fit$rotation), , drop = FALSE])
    ) |>
    as.data.frame() |>
    rownames_to_column("patient_id") |>
    select(patient_id, PC1)

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

  list(
    PCA = standardize_scores(pca_training, pca_test, "PC1"),
    pca_fit = pca_fit
  )
}
