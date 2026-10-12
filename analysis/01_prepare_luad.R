source("analysis/00_setup.R")

fetch <- function(url, path) {
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
    tmp <- paste0(path, '.partial')
    on.exit(unlink(tmp))
    utils::download.file(url, tmp, mode='wb', quiet=TRUE)
    if (file.info(tmp)$size == 0) stop('Empty download: ', url)
    if (!file.rename(tmp,path)) stop('Cannot save download: ',path)
  }
}
prepare_luad <- function() {
  prefix <- 'https://linkedomics.org/data_download/TCGA-LUAD/'
  fetch(paste0(prefix,'Human__TCGA_LUAD__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi'), 'data/raw/clinical.tsv')
  fetch(paste0(prefix,'Human__TCGA_LUAD__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz'), 'data/raw/expression.tsv.gz')
  raw <- readr::read_tsv('data/raw/clinical.tsv', show_col_types=FALSE)
  clinical <- data.frame(patient_id=names(raw)[-1],
    OS_time=as.numeric(unlist(raw[match('overall_survival',raw[[1]]),-1])),
    OS_event=as.numeric(unlist(raw[match('status',raw[[1]]),-1])))
  # Preserve original patient column names explicitly; tibble row names are not IDs.
  keep <- with(clinical,is.finite(OS_time) & OS_time>0 & OS_event %in% c(0,1))
  exclusion <- c(clinical_total=nrow(clinical),clinical_invalid=sum(!keep))
  clinical <- clinical[keep,]
  raw <- readr::read_tsv('data/raw/expression.tsv.gz',show_col_types=FALSE)
  genes <- raw[[1]]
  expression <- as.matrix(raw[,-1]); storage.mode(expression) <- 'double'
  # Collapse duplicate symbols by averaging their log2 expression, if present.
  if (anyDuplicated(genes)) {
    sums <- rowsum(expression,genes,reorder=FALSE)
    expression <- sums/as.numeric(table(factor(genes,levels=rownames(sums))))
  } else rownames(expression) <- genes
  patients <- clinical$patient_id[clinical$patient_id %in% colnames(expression)]
  exclusion <- c(exclusion,clinical_without_expression=sum(!clinical$patient_id %in% patients),matched=length(patients))
  clinical <- clinical[match(patients,clinical$patient_id),]
  expression <- expression[,patients,drop=FALSE]
  stopifnot(!anyDuplicated(clinical$patient_id), identical(colnames(expression),clinical$patient_id))
  dir.create('data/processed',recursive=TRUE,showWarnings=FALSE)
  saveRDS(list(clinical=clinical,expression=expression,exclusion=exclusion), 'data/processed/luad_analysis_input.rds')
  print(exclusion); print(table(clinical$OS_event))
  list(clinical=clinical,expression=expression,exclusion=exclusion)
}

prepare_luad_stage <- function(directory=stage_dir) save_stage(prepare_luad(), "input", directory)
if (sys.nframe()==0) prepare_luad_stage()
