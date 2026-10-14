library(tidyverse)
library(pathwayPCA)

# Load function
source("pathway_score_functions/03_pathwayPCA_function.R")


# Convert patient-by-gene LUAD data into the format your PCA helper expects.
run_pathway_pca <- function(gene_expression,
                            pathway_data) {

  # Identify the gene-expression columns
  gene_columns <- setdiff(
    names(gene_expression),
    c("patient_id", "OS_time", "OS_event")
  )

  # Put genes along rows and patients along columns
  expression_matrix <- as.matrix(
    gene_expression[, gene_columns, drop = FALSE]
  )

  rownames(expression_matrix) <- gene_expression$patient_id
  expression_matrix <- t(expression_matrix)

  # Exclude constant genes and genes with missing expression
  keep <- apply(
    expression_matrix,
    1,
    function(x) {
      all(is.finite(x)) && sd(x) > 0
    }
  )

  expression_data <- as_tibble(
    expression_matrix[keep, , drop = FALSE],
    rownames = "gene"
  )

  # Calculate PCA scores for each pathway
  results <- lapply(
    seq_along(pathway_data$pathways),
    function(i) {

      genes <- intersect(
        pathway_data$pathways[[i]],
        expression_data$gene
      )

      # Skip pathways with fewer than two usable genes
      if (length(genes) < 2) {
        return(NULL)
      }

      scores <- pathway_pca(
        expression_data,
        genes
      )

      scores$pathway <- pathway_data$TERMS[i]

      scores
    }
  )

  # Combine scores from all usable pathways
  results <- bind_rows(results)

  if (nrow(results) == 0) {
    stop("No pathways have at least two usable genes")
  }

  results
}


# Load LUAD data
luaddata <- readRDS(
  "data/processed/tcga_luad_processed.rds"
)


# Load KEGG pathways
kegg_pathways <- read_gmt(
  "data/pathways/kegg_medicus.gmt",
  description = TRUE
)


# Run PathwayPCA
kegg_results <- run_pathway_pca(
  gene_expression = luaddata,
  pathway_data = kegg_pathways
)


# Check results
dim(kegg_results)

head(kegg_results)

kegg_results