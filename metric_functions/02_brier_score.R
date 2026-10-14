library(survival)
library(purrr)

brier_score <- function(actual_time,
                        actual_event,
                        predicted_survival,
                        time_point) {
  
  # Estimate the censoring distribution
  censoring_fit <- survfit(
    Surv(actual_time, 1 - actual_event) ~ 1
  )
  
  # Function for obtaining G(t)
  get_G <- function(t,
                    left = FALSE) {

    # Use G(t-) for an event and G(t) at the evaluation time
    indices <- which(
      if (left) {
        censoring_fit$time < t
      } else {
        censoring_fit$time <= t
      }
    )

    if (length(indices) == 0) {
      return(1)
    }

    tail(
      censoring_fit$surv[indices],
      1
    )
  }

  # Calculate each participant's contribution
  contributions <- map_dbl(
    seq_along(actual_time),
    function(i) {
      
      # Event occurred by the evaluation time
      if (actual_time[i] <= time_point &&
          actual_event[i] == 1) {
        
        G_i <- get_G(actual_time[i], left = TRUE)
        
        return(
          predicted_survival[i]^2 / G_i
        )
      }
      
      # Participant survived beyond the evaluation time
      if (actual_time[i] > time_point) {
        
        G_t <- get_G(time_point)
        
        return(
          (1 - predicted_survival[i])^2 / G_t
        )
      }
      
      # Censored before the evaluation time
      return(0)
    }
  )
  
  # Average across participants
  score <- mean(contributions)
  
  return(score)
}
