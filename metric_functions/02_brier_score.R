# Estimate G from training outcomes; left limits are used at event times.
censoring_survival <- function(training_time, training_event) {
  fit <- survival::survfit(survival::Surv(training_time, 1-training_event) ~ 1)
  function(t, left=FALSE) vapply(t, function(x) {
    i <- which(if (left) fit$time < x else fit$time <= x)
    if (length(i)) tail(fit$surv[i], 1) else 1
  }, numeric(1))
}
brier_score <- function(actual_time, actual_event, predicted_survival, time_point,
                        training_time=actual_time, training_event=actual_event, G=NULL) {
  stopifnot(length(actual_time)==length(actual_event), length(predicted_survival)==length(actual_time),
            all(is.finite(predicted_survival)), all(predicted_survival>=0 & predicted_survival<=1))
  if (is.null(G)) G <- censoring_survival(training_time, training_event)
  event <- actual_time<=time_point & actual_event==1
  alive <- actual_time>time_point
  g_event <- G(actual_time[event], left=TRUE); g_t <- G(time_point)
  if (any(g_event<=0) || g_t<=0) stop("Censoring support exhausted")
  score <- numeric(length(actual_time))
  score[event] <- predicted_survival[event]^2/g_event
  score[alive] <- (1-predicted_survival[alive])^2/g_t
  mean(score)
}
