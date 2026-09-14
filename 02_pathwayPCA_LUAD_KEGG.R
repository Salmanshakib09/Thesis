library(pathwayPCA)
library(tidyverse)
# Read processed TCGA-LUAD data
luad_data <- readRDS(
  "data/processed/tcga_luad_processed.rds"
)

# Read KEGG pathways
kegg_pathways <- read_gmt(
  "data/pathways/kegg_medicus.gmt",
  description=TRUE
  )
luad_assay <- luad_data|>
select(-OS_time,-OS_event)
luad_survival= luad_data|>
  select(patient_id,OS_time,OS_event)

luad_omics <- CreateOmics(
  assayData_df = luad_assay,
  pathwayCollection_ls = kegg_pathways,
  response = luad_survival,
  respType = "surv"
)
class(luad_omics)
luad_omics
luad_aespc <- AESPCA_pVals(
  object = luad_omics,
  numPCs = 1,
  numReps = 0,
  parallel = TRUE,
  numCores = 2,
  adjustpValues = TRUE,
  adjustment = "BH"
)
top_paths= getPathpVals(luad_aespc, 
                        numPaths = 20
                        )
top_path <- top_paths$terms[1]

top_path
top_path_results <- getPathPCLs(
  luad_aespc,
  top_path
)
top_path_results