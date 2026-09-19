library(tidyverse)
library(pathwayPCA)

# Load function
source("03_pathwayPCA_function.R")


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