library(survival)
library(survAUC)

uno_cindex <- function(actual_time,
                       actual_event,
                       predicted_risk,
                       tau,
                       training_time = actual_time,
                       training_event = actual_event) {
  
  # Create the survival outcome
  actual_survival <- Surv(
    time = actual_time,
    event = actual_event
  )
  
  # Calculate Uno's C-index
  cindex <- UnoC(
    Surv.rsp = Surv(training_time, training_event),
    Surv.rsp.new = actual_survival,
    lpnew = predicted_risk,
    time = tau
  )
  
  return(as.numeric(cindex))
}
