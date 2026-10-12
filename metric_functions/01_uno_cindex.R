# Optional training outcomes supply the censoring reference for held-out tests.
uno_cindex <- function(actual_time, actual_event, predicted_risk, tau,
                       training_time=actual_time, training_event=actual_event) {
  as.numeric(survAUC::UnoC(
    survival::Surv(training_time, training_event),
    survival::Surv(actual_time, actual_event), lpnew=predicted_risk, time=tau))
}
