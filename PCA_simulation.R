# Clear old objects from R
rm(list = ls())


# Read KEGG pathway file
kegg_lines <- readLines(
  "data/raw/pathways/c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt"
)

length(kegg_lines)


# Split each pathway
kegg_split <- strsplit(kegg_lines, "\t")


# Get pathway names
pathway_names <- sapply(
  kegg_split,
  function(x) x[1]
)


# Get genes belonging to each pathway
kegg_pathways <- lapply(
  kegg_split,
  function(x) x[3:length(x)]
)

names(kegg_pathways) <- pathway_names


# Choose one KEGG pathway
selected_pathway <-
  "KEGG_MEDICUS_REFERENCE_CDC25_CELL_CYCLE_G2_M"


# Get genes from that pathway
selected_genes <- kegg_pathways[[selected_pathway]]

selected_genes

length(selected_genes)


# -------------------------
# SIMULATE EXPRESSION DATA
# -------------------------

set.seed(123)

expression_data <- matrix(
  rnorm(100 * length(selected_genes)),
  nrow = 100,
  ncol = length(selected_genes)
)


# Add gene names
colnames(expression_data) <- selected_genes


# Add patient names
rownames(expression_data) <- paste0(
  "Patient_",
  1:100
)


# Check simulated data
dim(expression_data)

expression_data[1:5, ]


# -------------------------
# PCA
# -------------------------

pca_result <- prcomp(
  expression_data,
  center = TRUE,
  scale. = TRUE
)


# PCA summary
summary(pca_result)


# -------------------------
# EIGENVALUES
# -------------------------

eigenvalues <- pca_result$sdev^2

eigenvalues


# Percentage variance explained
variance_explained <-
  eigenvalues / sum(eigenvalues) * 100

variance_explained


# Cumulative variance
cumulative_variance <-
  cumsum(variance_explained)

cumulative_variance


# -------------------------
# LOADINGS
# -------------------------

# Gene weights for PC1
pca_result$rotation[, 1]


# Gene weights for PC2
pca_result$rotation[, 2]


# -------------------------
# PATIENT PC SCORES
# -------------------------

pca_result$x[1:5, 1:2]


# -------------------------
# SCREE PLOT
# -------------------------

plot(
  variance_explained,
  type = "b",
  xlab = "Principal Component",
  ylab = "Percentage of Variance Explained",
  main = "Scree Plot"
)


# -------------------------
# PC1 VS PC2
# -------------------------

plot(
  pca_result$x[, 1],
  pca_result$x[, 2],
  xlab = "PC1",
  ylab = "PC2",
  main = "PCA of Simulated Pathway Expression Data"
)


# -------------------------
# SAVE RESULTS
# -------------------------

saveRDS(
  expression_data,
  "simulated_expression_data.rds"
)

saveRDS(
  pca_result,
  "pca_simulation_result.rds"
)