library(survival)
library(purrr)

brier_score <- function(actual_time,
                        actual_event,
                        predicted_survival,
                        time_point) {
  
  if (length(actual_time) == 0 ||
      length(actual_time) != length(actual_event) ||
      length(actual_time) != length(predicted_survival) ||
      any(!is.finite(actual_time)) ||
      anyNA(actual_event) ||
      any(!actual_event %in% c(0, 1)) ||
      any(!is.finite(predicted_survival)) ||
      any(predicted_survival < 0 | predicted_survival > 1) ||
      length(time_point) != 1 || !is.finite(time_point)) {
    stop("Provide aligned, finite outcomes and survival probabilities in [0, 1].")
  }

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
        if (G_i <= 0) {
          stop("Censoring survival is zero at an evaluated event time.")
        }
        
        return(
          predicted_survival[i]^2 / G_i
        )
      }
      
      # Participant survived beyond the evaluation time
      if (actual_time[i] > time_point) {
        
        G_t <- get_G(time_point)
        if (G_t <= 0) {
          stop("Censoring survival is zero at the evaluation horizon.")
        }
        
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
