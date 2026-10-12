source("metric_functions/02_brier_score.R")
integrated_brier_score <- function(actual_time, actual_event, predicted_survival, times,
                                   training_time=actual_time, training_event=actual_event) {
  if (length(times)<2 || any(diff(times)<=0)) stop("Use at least two increasing times")
  stopifnot(nrow(predicted_survival)==length(actual_time), ncol(predicted_survival)==length(times))
  G <- censoring_survival(training_time, training_event)
  scores <- vapply(seq_along(times), function(j) brier_score(actual_time, actual_event,
    predicted_survival[,j], times[j], G=G), numeric(1))
  sum(diff(times)*(head(scores,-1)+tail(scores,-1))/2)/diff(range(times))
}
