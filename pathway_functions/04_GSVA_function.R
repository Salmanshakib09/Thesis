# Keep all expression background genes. A vector gives one GSVA column;
# a named list gives one column per pathway. GSVA is cohort-dependent.
pathway_gsva <- function(expression_data, pathway_genes) {
  if (!"gene" %in% names(expression_data) || anyDuplicated(expression_data$gene)) stop("Supply unique gene identifiers")
  single <- !is.list(pathway_genes)
  sets <- if (single) list(GSVA=pathway_genes) else pathway_genes
  if (is.null(names(sets)) || anyDuplicated(names(sets))) stop("Use unique pathway names")
  sets <- lapply(sets, intersect, expression_data$gene)
  if (any(lengths(sets) < 2)) stop("Fewer than two pathway genes present")
  x <- as.matrix(expression_data[, setdiff(names(expression_data), "gene"), drop=FALSE])
  rownames(x) <- expression_data$gene
  if (any(!is.finite(x))) stop("GSVA expression must be finite")
  scores <- t(GSVA::gsva(GSVA::gsvaParam(x, sets, kcdf="Gaussian"), verbose=FALSE))
  if (any(!is.finite(scores))) stop("Non-finite GSVA scores: retain sufficient background genes")
  data.frame(patient_id=rownames(scores), scores, check.names=FALSE)
}
