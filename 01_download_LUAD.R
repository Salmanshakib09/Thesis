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
