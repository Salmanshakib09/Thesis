# Source the runner without executing it so all stages are available.
source("analysis/luad_reactome.R")
test_modular_workflow <- function() {
  set.seed(seed); n <- 160; p <- 400
  x <- matrix(rnorm(p*n),p,n,dimnames=list(paste0('g',1:p),paste0('patient',1:n)))
  signal <- colMeans(x[1:20,]); event_time <- rexp(n,exp(signal)); censor <- rexp(n,0.4)
  clinical <- data.frame(patient_id=colnames(x),OS_time=pmin(event_time,censor),OS_event=as.integer(event_time<=censor))
  sets <- setNames(lapply(1:12,function(i) paste0('g',((i-1)*25+1):((i-1)*25+25))),paste0('SYNTHETIC_',1:12))
  # A no-censoring hand calculation checks metric weights independently.
  check <- data.frame(OS_time=c(1,3),OS_event=c(1,1))
  stopifnot(abs(brier_score(check$OS_time,check$OS_event,c(0.2,0.8),2)-0.04)<1e-12)
  directory <- '/workspace/setup/luad-modular-self-test/stages'
  save_stage(list(clinical=clinical,expression=x,exclusion='SYNTHETIC VALIDATION ONLY'),'input',directory)
  eligible <- sort(names(sets)); set.seed(seed); selected <- sort(sample(eligible,10))
  save_stage(list(sets=sets[selected],version='synthetic'),'reactome',directory)
  run_stages(directory)
  cat('PASS: four configurations; finite metrics; survival probabilities; synthetic data only\n')

  # IPCW left-limit check: a censoring jump at time 1 must not weight an event at 1.
  G <- censoring_survival(c(1,1,3),c(1,0,1))
  stopifnot(G(1,left=TRUE)==1, abs(G(1)-2/3)<1e-12)
  # Test projection must not change the fitted PCA loadings.
  frame <- expression_frame(x[,1:30]); new <- expression_frame(x[,31:40])
  tr <- pathway_pca(frame,rownames(x)[1:25]); fit <- attr(tr,"pca_fit")
  te <- pathway_pca(new,rownames(x)[1:25],pca_fit=fit)
  stopifnot(identical(attr(te,"pca_fit"),fit),nrow(te)==10)
}
if (sys.nframe()==0) test_modular_workflow()
