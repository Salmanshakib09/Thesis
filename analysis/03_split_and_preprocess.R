source("analysis/00_setup.R")

split_patients <- function(clinical) {
  set.seed(seed)
  train <- sort(unlist(lapply(split(seq_len(nrow(clinical)),clinical$OS_event),function(i) sample(i,floor(0.7*length(i))))))
  list(train=train,test=setdiff(seq_len(nrow(clinical)),train))
}

split_stage <- function(directory=stage_dir) {
  input <- read_stage("input", directory); db <- read_stage("reactome", directory)
  split <- split_patients(input$clinical)
  training <- input$clinical[split$train,]; test <- input$clinical[split$test,]
  stopifnot(!length(intersect(training$patient_id,test$patient_id)))
  train <- input$expression[,split$train,drop=FALSE]; test_x <- input$expression[,split$test,drop=FALSE]
  med <- apply(train,1,median,na.rm=TRUE); valid <- is.finite(med)
  train <- train[valid,,drop=FALSE]; test_x <- test_x[valid,,drop=FALSE]; med <- med[valid]
  for (i in seq_len(nrow(train))) {
    train[i,!is.finite(train[i,])] <- med[i]
    test_x[i,!is.finite(test_x[i,])] <- med[i]
  }
  valid <- apply(train,1,sd)>0
  train <- train[valid,,drop=FALSE]; test_x <- test_x[valid,,drop=FALSE]
  sets <- lapply(db$sets,intersect,rownames(train))
  if (any(lengths(sets)<2)) stop("Insufficient genes after preprocessing")
  times <- seq(as.numeric(quantile(training$OS_time,0.1)),as.numeric(quantile(training$OS_time,0.6)),length.out=51)
  if (censoring_survival(training$OS_time,training$OS_event)(tail(times,1))<0.1) stop("Insufficient censoring support")
  save_stage(list(training=training,test=test,train_expression=train,test_expression=test_x,
    sets=sets,times=times,split=split,gene_medians=med,db_version=db$version),"split",directory)
}
if (sys.nframe()==0) split_stage()
