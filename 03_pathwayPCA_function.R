run_pathway_pca <- function(gene_expression, pathway_data) {
  
  # Separate gene expression data
  assay_data <- gene_expression |>
    dplyr::select(-OS_time, -OS_event)
  
  # Separate survival data
  survival_data <- gene_expression |>
    dplyr::select(patient_id, OS_time, OS_event)
  
  # Create PathwayPCA object
  omics_object <- pathwayPCA::CreateOmics(
    assayData_df = assay_data,
    pathwayCollection_ls = pathway_data,
    response = survival_data,
    respType = "surv"
  )
  
  # Run AES-PCA
  aespc_results <- pathwayPCA::AESPCA_pVals(
    object = omics_object,
    numPCs = 1,
    numReps = 0,
    parallel = TRUE,
    numCores = 2,
    adjustpValues = TRUE,
    adjustment = "BH"
  )
  
  # Extract results for ALL tested pathways
  results_table <- pathwayPCA::getPathpVals(
    aespc_results,
    numPaths = nrow(aespc_results$pVals_df)
  )
  
  return(results_table)
}