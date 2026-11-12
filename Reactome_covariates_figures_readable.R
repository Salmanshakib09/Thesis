# Run after Reactome_vector_covariates.R. Reads saved results; does not rerun GSVA.
# Each figure states its question, reference, evaluation sample, and limitations.

library(tidyverse)
library(survival)
source("metric_functions/01_uno_cindex.R")
source("metric_functions/02_brier_score.R")

input_folder <- "data/processed/reactome_vector_covariates"
figure_folder <- file.path(input_folder, "figures", "readable")
dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)

analysis <- readRDS(file.path(input_folder, "analysis.rds"))
training_data <- analysis$patient_split |>
  filter(split == "training")
test_data <- analysis$patient_split |>
  filter(split == "test")
pathway_models <- analysis$models
selected_pathways <- analysis$selected_pathways

stopifnot(
  length(pathway_models) == 10,
  !anyDuplicated(analysis$patient_split$patient_id),
  !any(training_data$patient_id %in% test_data$patient_id)
)

ridge_file <- file.path(input_folder, "ridge_analysis.rds")
ridge_result <- if (file.exists(ridge_file)) readRDS(ridge_file) else NULL

# Use the ridge evaluation grid if available; otherwise the same training rule.
evaluation_times <- if (!is.null(ridge_result)) {
  ridge_result$times
} else {
  seq(
    quantile(training_data$OS_time, 0.1),
    quantile(training_data$OS_time, 0.6),
    length.out = 51
  )
}

cohort_note <- paste0(
  "Training: n = ", nrow(training_data), ", deaths = ", sum(training_data$OS_event),
  "; test: n = ", nrow(test_data), ", deaths = ", sum(test_data$OS_event),
  ". Excluded for covariates: ", nrow(analysis$excluded_patients), ". One fixed split."
)
time_note <- sprintf(
  "Brier at %.1f days; IBS over %.1f to %.1f days.",
  last(evaluation_times), first(evaluation_times), last(evaluation_times)
)

pathway_label <- function(x) {
  x |>
    str_remove("^REACTOME_") |>
    str_replace_all("_", " ") |>
    str_wrap(width = 30)
}

# Explanatory notes belong in the QMD captions, not inside the images.
# Wrap titles and subtitles and use a consistent readable export layout.
save_figure <- function(plot, filename, width = 10, height = 8) {
  wrap_note <- function(note) {
    if (is.null(note)) return(NULL)
    str_split(note, "\\n")[[1]] |>
      str_wrap(width = floor(width * 9)) |>
      paste(collapse = "\n")
  }
  plot <- plot +
    labs(
      title = str_wrap(plot$labels$title, width = floor(width * 7)),
      subtitle = wrap_note(plot$labels$subtitle),
      caption = NULL
    ) +
    theme(
      text = element_text(size = 16),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 13),
      plot.title = element_text(size = 17, face = "bold"),
      plot.subtitle = element_text(size = 13, lineheight = 1.15,
                                   margin = margin(b = 12)),
      plot.caption = element_text(size = 10, hjust = 0, lineheight = 1.2,
                                  margin = margin(t = 14)),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.margin = margin(16, 20, 16, 16),
      panel.spacing = grid::unit(1, "lines"),
      legend.position = "bottom"
    )
  ggsave(file.path(figure_folder, paste0(filename, "_readable.png")),
         plot = plot, width = width, height = height, dpi = 300)
  ggsave(file.path(figure_folder, paste0(filename, "_readable.pdf")),
         plot = plot, width = width, height = height)
  invisible(plot)
}


# 1. Association: which score coefficients differ from a hazard ratio of one? ----

score_coefficients <- pathway_models |>
  imap(function(result, pathway) {
    estimates <- summary(result$model)$coefficients
    intervals <- summary(result$model)$conf.int
    tibble(
      pathway = pathway,
      term = rownames(estimates),
      hazard_ratio = intervals[, "exp(coef)"],
      lower_95 = intervals[, "lower .95"],
      upper_95 = intervals[, "upper .95"],
      p_value = estimates[, "Pr(>|z|)"]
    )
  }) |>
  list_rbind() |>
  filter(term %in% c("PC1", "GSVA", "ssGSEA")) |>
  mutate(
    adjusted_p = p.adjust(p_value, method = "BH"),
    term = factor(term, levels = c("PC1", "GSVA", "ssGSEA")),
    pathway_name = pathway_label(pathway)
  )

forest_plot <- score_coefficients |>
  ggplot(aes(hazard_ratio, pathway_name)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_segment(aes(x = lower_95, xend = upper_95, yend = pathway_name)) +
  geom_point(size = 2) +
  scale_x_log10() +
  facet_wrap(~ term, nrow = 1) +
  labs(
    title = "How are pathway scores associated with mortality hazard?",
    subtitle = paste0("Training associations only. n = ", nrow(training_data),
                      "; observed deaths = ", sum(training_data$OS_event)),
    x = "Adjusted hazard ratio per training SD (log scale)", y = NULL,
    caption = paste0(
      "Dashed reference: HR = 1. Individual 95% intervals; adjusted for the other two scores, age, sex and stage.\n",
      "BH adjustment across ", nrow(score_coefficients), " score tests: ",
      sum(score_coefficients$adjusted_p < 0.05), " pass 5% FDR. ",
      "PC1 sign is arbitrary. Associations do not demonstrate predictive improvement."
    )
  ) +
  theme_bw()

save_figure(forest_plot, "01_score_associations", height = 8)
score_coefficients |> write_csv(file.path(figure_folder, "score_associations.csv"))


# 2. Redundancy: do the three scores contain similar expression information? ----

correlations <- pathway_models |>
  imap(function(result, pathway) {
    scores <- result$training_scores
    tibble(
      pathway = pathway,
      comparison = c("PC1 vs GSVA", "PC1 vs ssGSEA", "GSVA vs ssGSEA"),
      correlation = c(
        cor(scores$PC1, scores$GSVA),
        cor(scores$PC1, scores$ssGSEA),
        cor(scores$GSVA, scores$ssGSEA)
      )
    )
  }) |>
  list_rbind() |>
  mutate(pathway_name = pathway_label(pathway))

correlation_plot <- correlations |>
  ggplot(aes(comparison, pathway_name, fill = correlation)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 4.5) +
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "firebrick",
    midpoint = 0, limits = c(-1, 1)
  ) +
  labs(
    title = "Do pathway-scoring methods provide overlapping information?",
    subtitle = paste0("Pearson correlations across ", nrow(training_data), " training patients"),
    x = NULL, y = NULL, fill = "Pearson r",
    caption = paste(
      "Reference: r = 0 means no linear correlation; values near either -1 or +1 indicate strong linear dependence.",
      "PC1 sign is arbitrary. No survival outcomes are used; this is not a prediction-accuracy comparison.",
      sep = "\n"
    )
  ) +
  theme_minimal()

save_figure(correlation_plot, "02_score_correlations", width = 10, height = 8)
correlations |> write_csv(file.path(figure_folder, "score_correlations.csv"))


# 3. Prediction: compare models against clinical-only Cox on identical patients ----

clinical_model <- training_data |>
  coxph(
    Surv(OS_time, OS_event) ~ age + sex + stage,
    data = _, na.action = na.fail, model = TRUE, x = TRUE
  )

# One prediction convention for clinical, pathway, and ridge Cox models.
cox_predictions <- function(model, newdata) {
  risk <- predict(model, newdata = newdata, type = "lp", reference = "sample") |>
    as.numeric()
  baseline <- basehaz(model, centered = TRUE)
  hazard <- evaluation_times |>
    map_dbl(function(t) {
      baseline |> filter(time <= t) |> pull(hazard) |> last(default = 0)
    })
  survival <- exp(-outer(exp(risk), hazard))
  stopifnot(length(risk) == nrow(newdata), all(is.finite(survival)))
  list(risk = risk, survival = survival)
}

evaluate_predictions <- function(outcomes, prediction, model_name, sample_name) {
  curve <- evaluation_times |>
    seq_along() |>
    map_dbl(function(j) {
      brier_score(outcomes$OS_time, outcomes$OS_event,
                  prediction$survival[, j], evaluation_times[j])
    })
  # Average over the stated interval, using the trapezoidal rule.
  ibs <- sum(diff(evaluation_times) * (head(curve, -1) + tail(curve, -1)) / 2) /
    diff(range(evaluation_times))
  list(
    metrics = tibble(
      model = model_name, sample = sample_name,
      n = nrow(outcomes), deaths = sum(outcomes$OS_event),
      uno_cindex = uno_cindex(
        outcomes$OS_time, outcomes$OS_event, prediction$risk,
        tau = last(evaluation_times), training_time = training_data$OS_time,
        training_event = training_data$OS_event
      ),
      brier_score = last(curve), integrated_brier_score = ibs
    ),
    curve = tibble(model = model_name, sample = sample_name,
                   time = evaluation_times, brier_score = curve)
  )
}

clinical_results <- list(
  evaluate_predictions(training_data, cox_predictions(clinical_model, training_data),
                       "Clinical-only Cox", "Training (apparent)"),
  evaluate_predictions(test_data, cox_predictions(clinical_model, test_data),
                       "Clinical-only Cox", "Held-out test")
)

pathway_results <- pathway_models |>
  imap(function(result, pathway) {
    stopifnot(
      identical(result$training_scores$patient_id, training_data$patient_id),
      identical(result$test_scores$patient_id, test_data$patient_id)
    )
    train <- training_data |> left_join(result$training_scores, by = "patient_id")
    test <- test_data |> left_join(result$test_scores, by = "patient_id")
    list(
      evaluate_predictions(training_data, cox_predictions(result$model, train),
                           pathway, "Training (apparent)"),
      evaluate_predictions(test_data, cox_predictions(result$model, test),
                           pathway, "Held-out test")
    )
  }) |>
  flatten()

performance <- c(clinical_results, pathway_results) |>
  map("metrics") |>
  list_rbind()

metric_labels <- c(
  uno_cindex = "Uno's C-index\n(higher is better)",
  brier_score = "Brier score\n(lower is better)",
  integrated_brier_score = "IBS\n(lower is better)"
)

performance_long <- performance |>
  pivot_longer(all_of(names(metric_labels)), names_to = "metric", values_to = "value") |>
  mutate(
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels),
    sample = factor(sample, levels = c("Training (apparent)", "Held-out test"))
  )

reference_values <- performance_long |>
  filter(model == "Clinical-only Cox")

pathway_performance_plot <- performance_long |>
  filter(model != "Clinical-only Cox") |>
  mutate(pathway_name = pathway_label(model)) |>
  ggplot(aes(value, pathway_name)) +
  geom_vline(data = reference_values, aes(xintercept = value),
             linetype = "dashed", colour = "firebrick", inherit.aes = FALSE) +
  geom_point(size = 2) +
  facet_grid(sample ~ metric, scales = "free_x") +
  labs(
    title = "Does each pathway model predict better than clinical information alone?",
    subtitle = paste(cohort_note, time_note, sep = "\n"),
    x = "Metric value", y = NULL,
    caption = paste(
      "Red dashed line: clinical-only Cox (age, sex, stage), fitted on the same training patients.",
      "Each point: one pathway model with PC1 + GSVA + ssGSEA + age + sex + stage.",
      "Training values are apparent, not cross-validated. Test results are point estimates from one split; no uncertainty intervals.",
      sep = "\n"
    )
  ) +
  theme_bw()

pathway_performance_plot <- save_figure(
  pathway_performance_plot, "03_pathway_vs_clinical", width = 15, height = 14
)


# Separate report pages avoid shrinking twenty pathway rows onto one page.
for (sample_name in levels(performance_long$sample)) {
  page_data <- performance_long |>
    filter(model != "Clinical-only Cox", sample == sample_name) |>
    mutate(pathway_name = pathway_label(model))
  page_reference <- reference_values |> filter(sample == sample_name)
  report_plot <- ggplot(page_data, aes(value, pathway_name)) +
    geom_vline(data = page_reference, aes(xintercept = value),
               linetype = "dashed", colour = "firebrick", inherit.aes = FALSE) +
    geom_point(size = 2) +
    facet_wrap(~ metric, scales = "free_x", nrow = 1) +
    scale_x_continuous(n.breaks = 3, labels = scales::label_number(accuracy = 0.001)) +
    labs(title = paste("Pathway versus clinical-only prediction:", sample_name),
         subtitle = time_note, x = "Metric value", y = NULL) +
    theme_bw()
  suffix <- if (sample_name == "Training (apparent)") "training" else "test"
  save_figure(report_plot, paste0("03_pathway_vs_clinical_", suffix),
              width = 10, height = 8)
}


# 4. Optional final-ridge figures, only when the current analysis saved that stage ----

if (!is.null(ridge_result)) {
  stopifnot(
    identical(ridge_result$training_data$patient_id, training_data$patient_id),
    identical(ridge_result$test_data$patient_id, test_data$patient_id),
    ncol(ridge_result$stack$oof) == length(pathway_models),
    identical(colnames(ridge_result$stack$oof), names(pathway_models))
  )

  ridge_results <- list(
    evaluate_predictions(
      training_data,
      cox_predictions(ridge_result$model, ridge_result$stack$training),
      "Ten-pathway stacked ridge", "Training (apparent)"
    ),
    evaluate_predictions(
      test_data,
      cox_predictions(ridge_result$model, ridge_result$stack$validation),
      "Ten-pathway stacked ridge", "Held-out test"
    )
  )

  final_results <- c(clinical_results, ridge_results)
  final_performance <- final_results |> map("metrics") |> list_rbind() |>
    mutate(sample = factor(sample, levels = c("Training (apparent)", "Held-out test")))
  final_curves <- final_results |> map("curve") |> list_rbind() |>
    mutate(sample = factor(sample, levels = c("Training (apparent)", "Held-out test")))

  final_plot <- final_performance |>
    pivot_longer(all_of(names(metric_labels)), names_to = "metric", values_to = "value") |>
    mutate(
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels),
    sample = factor(sample, levels = c("Training (apparent)", "Held-out test"))
  ) |>
    ggplot(aes(model, value, colour = sample, group = sample)) +
    geom_point(position = position_dodge(width = 0.4), size = 3) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(
      title = "Does the final stacked model improve on clinical-only prediction?",
      subtitle = paste(cohort_note, time_note, sep = "\n"),
      x = NULL, y = "Metric value", colour = "Evaluation sample",
      caption = paste(
        "Reference model: clinical-only Cox on identical patients. Clinical information also enters all first-stage pathway models.",
        "Ridge training metrics fit the meta-model on OOF features; they are still apparent, not OOF evaluation of the entire stack.",
        "One test split; no confidence intervals. This comparison does not establish improvement over gene-level or random-pathway models.",
        sep = "\n"
      )
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
    scale_x_discrete(labels = c("Clinical-only Cox" = "Clinical-only\nCox",
                                "Ten-pathway stacked ridge" = "Stacked\nridge"))

  save_figure(final_plot, "04_stacked_vs_clinical", width = 10, height = 8)

  brier_plot <- final_curves |>
    ggplot(aes(time, brier_score, colour = model, linetype = model)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ sample) +
    labs(
      title = "How does survival-probability prediction error change over follow-up?",
      subtitle = paste(cohort_note, "Lower Brier scores are better.", sep = "\n"),
      x = "Days since the source follow-up origin", y = "IPCW Brier score",
      colour = "Model", linetype = "Model",
      caption = paste(
        "Clinical-only Cox is the reference; both models use the same evaluation grid and censoring-weight convention.",
        "Training is apparent performance. A lower curve on this split is not a significance test. GSVA remains cohort-dependent.",
        sep = "\n"
      )
    ) +
    theme_bw()

  save_figure(brier_plot, "05_brier_curves", width = 10, height = 8)
  final_performance |> write_csv(file.path(figure_folder, "stacked_vs_clinical_metrics.csv"))
  final_curves |> write_csv(file.path(figure_folder, "stacked_vs_clinical_brier_curves.csv"))
  performance <- bind_rows(performance, final_performance |> filter(model != "Clinical-only Cox"))
} else {
  message(
    "Figures 1–3 saved. Ridge figures 4–5 were NOT generated: ridge_analysis.rds is absent. ",
    "Run the current Reactome_vector_covariates.R through its final ridge stage, then rerun this plotting script."
  )
}

performance |> write_csv(file.path(figure_folder, "model_metrics.csv"))
selected_pathways |> write_csv(file.path(figure_folder, "pathway_gene_counts.csv"))

tibble(
  figure = c("01_score_associations", "02_score_correlations", "03_pathway_vs_clinical",
             "04_stacked_vs_clinical", "05_brier_curves"),
  question = c(
    "Which score coefficients show conditional associations with mortality hazard?",
    "Do the three scores provide redundant expression information?",
    "Do individual pathway models improve on clinical-only prediction?",
    "Does the full stack improve on clinical-only prediction?",
    "How does probability prediction error vary with follow-up?"
  ),
  reference = c("HR = 1", "Pearson r = 0", "Clinical-only Cox", "Clinical-only Cox", "Clinical-only Cox"),
  generated = c(TRUE, TRUE, TRUE, !is.null(ridge_result), !is.null(ridge_result))
) |>
  write_csv(file.path(figure_folder, "figure_questions.csv"))

message("Figures and supporting tables saved under: ", figure_folder)
if (interactive()) pathway_performance_plot |> print()
