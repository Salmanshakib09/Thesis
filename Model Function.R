library(tidyverse)
library(survival)


fit_survival_model <- function(survival_data,
                               pathway_scores,
                               model = c("cox", "aft")) {
  
  model <- match.arg(model)
  
  # Identify the pathway-score column
  score_name <- pathway_scores |>
    select(-patient_id) |>
    names()
  
  if (length(score_name) != 1) {
    stop("pathway_scores must contain patient_id and one pathway-score column.")
  }
  
  # Rename the pathway score to a common name
  pathway_scores_clean <- pathway_scores |>
    rename(pathway_score = all_of(score_name))
  
  # Join survival data with pathway scores
  model_data <- survival_data |>
    inner_join(
      pathway_scores_clean,
      by = "patient_id"
    )
  
  # Fit Cox PH model
  if (model == "cox") {
    
    fitted_model <- model_data |>
      coxph(
        formula = Surv(OS_time, OS_event) ~ pathway_score,
        data = _
      )
    
  }
  
  # Fit Weibull AFT model
  if (model == "aft") {
    
    fitted_model <- model_data |>
      survreg(
        formula = Surv(OS_time, OS_event) ~ pathway_score,
        data = _,
        dist = "weibull"
      )
    
  }
  
  # Create prediction function
  predict_function <- function(new_pathway_scores) {
    
    new_data <- new_pathway_scores |>
      rename(
        pathway_score = all_of(score_name)
      )
    
    predictions <- fitted_model |>
      predict(
        newdata = new_data,
        type = "lp"
      )
    
    return(predictions)
  }
  
  return(predict_function)
}