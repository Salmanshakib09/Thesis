library(tidyverse)
library(survival)


fit_survival_model <- function(survival_data,
                               pathway_scores,
                               model = c("cox", "aft")) {
  
  model <- match.arg(model)
  
  # Identify all pathway-score columns
  score_names <- pathway_scores |>
    select(-patient_id) |>
    names()

  if (length(score_names) < 1) {
    stop("Provide patient_id and at least one pathway-score column.")
  }

  if (anyDuplicated(survival_data$patient_id) ||
      anyDuplicated(pathway_scores$patient_id)) {
    stop("Patient IDs must be unique.")
  }

  if (!all(survival_data$patient_id %in% pathway_scores$patient_id)) {
    stop("Some training patients have no pathway scores.")
  }

  # Join all score columns to the survival outcomes
  model_data <- survival_data |>
    inner_join(
      pathway_scores,
      by = "patient_id"
    )

  # Each pathway gets its own coefficient in the SAME model
  model_formula <- score_names |>
    paste0("`", x = _, "`") |>
    reformulate(
      response = "Surv(OS_time, OS_event)"
    )

  # Fit Cox PH model
  if (model == "cox") {
    
    fitted_model <- model_data |>
      coxph(
        formula = model_formula,
        data = _
      )
    
  }
  
  # Fit Weibull AFT model
  if (model == "aft") {
    
    fitted_model <- model_data |>
      survreg(
        formula = model_formula,
        data = _,
        dist = "weibull"
      )
    
  }
  
  # Create prediction function
  predict_function <- function(new_pathway_scores) {
    
    new_data <- new_pathway_scores |>
      select(all_of(score_names))

    predictions <- fitted_model |>
      predict(
        newdata = new_data,
        type = "lp"
      )
    
    return(predictions)
  }
  
  return(predict_function)
}