library(survival)
library(survAUC)

uno_cindex <- function(actual_time,
                       actual_event,
                       predicted_risk,
                       tau) {
  
  # Create the survival outcome
  actual_survival <- Surv(
    time = actual_time,
    event = actual_event
  )
  
  # Calculate Uno's C-index
  cindex <- UnoC(
    Surv.rsp = actual_survival,
    Surv.rsp.new = actual_survival,
    lpnew = predicted_risk,
    time = tau
  )
  
  return(as.numeric(cindex))
}