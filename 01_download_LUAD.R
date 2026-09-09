library(readr)

clinical_url <- "https://linkedomics.org/data_download/TCGA-LUAD/Human__TCGA_LUAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi"

clinical_raw <- read_tsv(clinical_url)

dim(clinical_raw)
head(clinical_raw)
clinical_raw$attrib_name
clinical_raw[clinical_raw$attrib_name %in%
               c("overall_survival", "status", "overallsurvival"),
             1:8]
# Patient IDs
patient_id <- colnames(clinical_raw)[-1]

# Overall survival time
OS_time <- as.numeric(
  unlist(
    clinical_raw[
      clinical_raw$attrib_name == "overall_survival",
      -1
    ]
  )
)

# Overall survival event
OS_event <- as.numeric(
  unlist(
    clinical_raw[
      clinical_raw$attrib_name == "status",
      -1
    ]
  )
)

# Create processed clinical dataset
clinical_processed <- data.frame(
  patient_id = patient_id,
  OS_time = OS_time,
  OS_event = OS_event
)
dim(clinical_processed)
head(clinical_processed)
colSums(is.na(clinical_processed))

table(clinical_processed$OS_event, useNA = "ifany")

summary(clinical_processed$OS_time)
clinical_clean <- clinical_processed[
  !is.na(clinical_processed$OS_time) &
    !is.na(clinical_processed$OS_event),
]
dim(clinical_clean)

table(clinical_clean$OS_event)
# RNA-seq gene expression data
expression_url <- "https://linkedomics.org/data_download/TCGA-LUAD/Human__TCGA_LUAD__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz"

expression_raw <- read_tsv(expression_url)

dim(expression_raw)
head(expression_raw[, 1:6])
# Patient IDs in the expression data
expression_patient_ids <- colnames(expression_raw)[-1]

# Patients available in both clinical and expression data
common_patients <- intersect(
  clinical_clean$patient_id,
  expression_patient_ids
)

length(common_patients)
# Keep only the 492 common patients in the clinical data
clinical_matched <- clinical_clean[
  match(common_patients, clinical_clean$patient_id),
]

# Keep only the same 492 patients in the expression data
expression_matched <- expression_raw[
  ,
  c("attrib_name", common_patients)
]
dim(clinical_matched)
dim(expression_matched)
identical(
  clinical_matched$patient_id,
  colnames(expression_matched)[-1]
)
# Any duplicated gene names?
sum(duplicated(expression_matched$attrib_name))

# Any missing expression values?
sum(is.na(expression_matched[, -1]))
# Convert expression data to a numeric matrix
expr_matrix <- as.matrix(expression_matched[, -1])

# Put gene names on the rows
rownames(expr_matrix) <- expression_matched$attrib_name

# Transpose so rows = patients and columns = genes
expr_patient_gene <- t(expr_matrix)
luad_processed <- data.frame(
  clinical_matched,
  expr_patient_gene,
  check.names = FALSE
)
dim(luad_processed)
head(luad_processed[, 1:8])
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(
  luad_processed,
  file.path("data", "processed", "tcga_luad_processed.rds")
)
test_data <- readRDS("data/processed/tcga_luad_processed.rds")

dim(test_data)
file.exists("data/processed/tcga_luad_processed.rds")