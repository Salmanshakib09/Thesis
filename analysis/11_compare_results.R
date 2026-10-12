source("analysis/00_setup.R")
compare_stage <- function(directory=stage_dir,outdir=dirname(directory)) {
  input <- read_stage("input",directory); x <- read_stage("split",directory)
  db <- list(version=x$db_version); clinical <- input$clinical
  training <- x$training; test <- x$test; times <- x$times; horizon <- tail(times,1)
  split <- x$split; selected <- names(x$sets)
  pca <- read_stage("PCA_scores",directory); gsva <- read_stage("GSVA_scores",directory)
  scored <- list(scores=list(PCA=pca$scores,GSVA=gsva$scores),pca_fits=pca$pca_fits,sets=x$sets,gene_medians=x$gene_medians)
  ci <- read_stage("uno_cindex",directory); bs <- read_stage("brier_score",directory); ibs <- read_stage("integrated_brier_score",directory)
  all_metrics <- curves <- predictions <- models <- list()
  for (key in names(ci)) {
    parts <- strsplit(key,"_")[[1]]; pred <- read_stage(key,directory)
    all_metrics[[key]] <- data.frame(scoring=parts[1],model=parts[2],n_pathways=length(selected),train_n=nrow(training),test_n=nrow(test),test_events=sum(test$OS_event),horizon=horizon,uno_cindex=ci[[key]],brier_score=tail(bs[[key]],1),integrated_brier_score=ibs[[key]])
    curves[[key]] <- data.frame(configuration=key,time=times,brier_score=bs[[key]])
    predictions[[key]] <- pred[c("patient_id","risk","survival")]; models[[key]] <- pred$model
  }
  dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
  metrics <- do.call(rbind,all_metrics); rownames(metrics) <- NULL
  stopifnot(nrow(metrics)==4,all(is.finite(as.matrix(metrics[,c('uno_cindex','brier_score','integrated_brier_score')]))),all(metrics$uno_cindex>=0 & metrics$uno_cindex<=1))
  write.csv(metrics,file.path(outdir,'performance.csv'),row.names=FALSE)
  write.csv(do.call(rbind,curves),file.path(outdir,'brier_curves.csv'),row.names=FALSE)
  write.csv(data.frame(feature=paste0('p',1:10),pathway=selected,measured_genes=lengths(scored$sets)),file.path(outdir,'selected_pathways.csv'),row.names=FALSE)
  write.csv(data.frame(clinical,split=ifelse(seq_len(nrow(clinical)) %in% split$train,'train','test')),file.path(outdir,'patient_split.csv'),row.names=FALSE)
  saveRDS(list(seed=seed,db_version=db$version,models=models,scored=scored,predictions=predictions,times=times,exclusion=input$exclusion),file.path(outdir,'analysis.rds'))
  capture.output(sessionInfo(),file=file.path(outdir,'session-info.txt'))
  capture.output(print(input$exclusion),file=file.path(outdir,'cohort-flow.txt'))
  for (key in c('PCA_Cox','GSVA_Cox')) {
    diagnostic <- survival::cox.zph(models[[key]])
    write.csv(diagnostic$table,file.path(outdir,paste0(key,'_PH_diagnostics.csv')))
  }
  if (!identical(db$version,'synthetic')) {
    inputs <- c('data/raw/clinical.tsv','data/raw/expression.tsv.gz','data/raw/reactome-msigdbr.rds', list.files('analysis',pattern='[.]R$',full.names=TRUE),list.files('pathway_functions',pattern='[.]R$',full.names=TRUE),list.files('metric_functions',pattern='[.]R$',full.names=TRUE),'Model Function.R')
    writeLines(system2('sha256sum',shQuote(inputs),stdout=TRUE),file.path(outdir,'input-checksums.sha256'))
  }
  grDevices::pdf(file.path(outdir,'brier_curves.pdf'),width=7,height=5)
  cols <- c('#0072B2','#D55E00','#009E73','#CC79A7')
  plot(range(times),range(unlist(lapply(curves,function(x) x$brier_score))),type='n',xlab='Survival time (source units)',ylab='IPCW Brier score',main='TCGA-LUAD: held-out prediction error')
  for(i in seq_along(curves)) lines(curves[[i]]$time,curves[[i]]$brier_score,col=cols[i],lwd=2)
  legend('topleft',names(curves),col=cols,lwd=2,bty='n'); dev.off()
  print(metrics); invisible(metrics)
}
if (sys.nframe()==0) compare_stage()
