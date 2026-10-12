# Accept one or more pathway-score columns, with explicit patient IDs.
# The returned function preserves the original one-argument LP prediction API.
fit_survival_model <- function(survival_data, pathway_scores, model=c("cox", "aft")) {
  model <- match.arg(model)
  score_names <- setdiff(names(pathway_scores), "patient_id")
  if (!length(score_names)) stop("At least one pathway score is required")
  if (anyDuplicated(survival_data$patient_id) || anyDuplicated(pathway_scores$patient_id)) stop("Duplicate patient IDs")
  i <- match(survival_data$patient_id, pathway_scores$patient_id)
  if (anyNA(i)) stop("Missing pathway scores for training patients")
  features <- paste0("p", seq_along(score_names))
  scores <- pathway_scores[i, score_names, drop=FALSE]; names(scores) <- features
  dat <- data.frame(survival_data, scores)
  formula <- as.formula(paste("survival::Surv(OS_time, OS_event) ~", paste(features, collapse="+")))
  fit <- if (model == "cox") survival::coxph(formula, data=dat, x=TRUE, y=TRUE) else
    survival::survreg(formula, data=dat, dist="weibull", x=TRUE, y=TRUE)
  if (any(!is.finite(coef(fit)))) stop("Non-finite model coefficients")
  prediction <- function(new_pathway_scores, times=NULL, type=c("lp", "risk", "survival")) {
    type <- match.arg(type)
    dat <- new_pathway_scores[, score_names, drop=FALSE]; names(dat) <- features
    lp <- if (model == "cox") as.numeric(predict(fit, newdata=dat, type="lp", reference=if (type=="lp") "sample" else "zero")) else
      as.numeric(predict(fit, newdata=dat, type="lp"))
    if (type == "lp") return(lp)
    if (type == "risk") return(if (model == "cox") lp else -lp)
    if (is.null(times) || any(!is.finite(times)) || any(times <= 0)) stop("Supply positive prediction times")
    if (model == "cox") {
      bh <- survival::basehaz(fit, centered=FALSE)
      h0 <- vapply(times, function(t) { i <- which(bh$time <= t); if (length(i)) tail(bh$hazard[i], 1) else 0 }, numeric(1))
      exp(-outer(exp(lp), h0))
    } else exp(-exp(outer(-lp, log(times), "+") / fit$scale))
  }
  attr(prediction, "model") <- fit
  prediction
}
