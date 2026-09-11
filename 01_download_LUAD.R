library(readr)
library(dplyr)

clinical_url <- "https://linkedomics.org/data_download/TCGA-LUAD/Human__TCGA_LUAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi"
clinical_raw <- read_tsv(clinical_url)

# Transpose the dataset to get partients along rows.
clinical_t = t(clinical_raw[, -1])
colnames(clinical_t) = clinical_raw[[1]]
clinical_t = as_tibble(clinical_t)
clinical_t$patient_id = rownames(clinical_t)

# Subset onto required variables.
clinical_processed = clinical_t |>
  mutate(OS_time = as.numeric(overall_survival),
         OS_event = as.numeric(status)) |>
  select(patient_id, OS_time, OS_event)

clinical_clean = clinical_processed |>
  filter(!is.na(OS_time),
         !is.na(OS_event))

# RNA-seq gene expression data
expression_url <- "https://linkedomics.org/data_download/TCGA-LUAD/Human__TCGA_LUAD__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz"
expression_raw <- read_tsv(expression_url)

# Transpose the dataset to get partients along rows.
expression_t = t(expression_raw[, -1])
colnames(expression_t) = expression_raw[[1]]
expression_t = as_tibble(expression_t)
expression_t$patient_id = rownames(expression_t)

# Join clinical data with gene expression data.
luad_processed = inner_join(clinical_clean, expression_t, by = "patient_id")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(
  luad_processed,
  file.path("data", "processed", "tcga_luad_processed.rds")
)