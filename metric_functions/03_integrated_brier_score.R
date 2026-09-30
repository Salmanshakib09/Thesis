library(purrr)

source("metric_functions/02_brier_score.R")

integrated_brier_score <- function(actual_time,
                                   actual_event,
                                   predicted_survival,
                                   times) {
  
  # Calculate Brier score at each evaluation time
  brier_scores <- map_dbl(
    seq_along(times),
    function(j) {
      
      brier_score(
        actual_time = actual_time,
        actual_event = actual_event,
        predicted_survival = predicted_survival[, j],
        time_point = times[j]
      )
    }
  )
  
  # Integrate the Brier scores using the trapezoidal rule
  area <- sum(
    diff(times) *
      (head(brier_scores, -1) +
         tail(brier_scores, -1)) / 2
  )
  
  # Calculate IBS
  ibs <- area / (max(times) - min(times))
  
  return(ibs)
}