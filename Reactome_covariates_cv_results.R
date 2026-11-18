# Run after Reactome_covariates_cross_validation.R; no model refitting here.
# Writes four figures, two tables, and an editable Results draft.
library(tidyverse)
library(survival)

input_folder <- "data/processed/reactome_covariates_cv"
figure_folder <- "figures/reactome_cv"
report_folder <- "results/reactome_cv"
dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(report_folder, recursive = TRUE, showWarnings = FALSE)
a <- readRDS(file.path(input_folder, "cv_analysis.rds"))
metrics <- a$metrics
stopifnot(a$settings$k == 5, a$settings$inner_folds == 3,
          nrow(metrics) == 60, "Stacked ridge" %in% metrics$model)
fixed_pathways <- a$settings$pathways
metric_names <- c("uno_cindex", "brier_score", "integrated_brier_score")
metric_labels <- c("Uno's C-index\n(higher is better)",
                   "Brier score\n(lower is better)", "IBS\n(lower is better)")
ids <- set_names(sprintf("P%02d", seq_along(fixed_pathways)), fixed_pathways)
pathway_names <- fixed_pathways |> str_remove("^REACTOME_") |> str_replace_all("_", " ")
labels <- set_names(paste(ids, pathway_names, sep = ": "), fixed_pathways)
labels <- c(labels, "Clinical-only Cox" = "Clinical-only Cox", "Stacked ridge" = "Stacked ridge")
model_order <- c(fixed_pathways, "Stacked ridge", "Clinical-only Cox")
model_label <- function(x) str_wrap(unname(labels[as.character(x)]), 32)

long <- metrics |> pivot_longer(all_of(metric_names), names_to = "metric", values_to = "value")
summary_metrics <- long |> group_by(model, metric) |>
  summarise(mean = mean(value), sd = sd(value), n_folds = n(),
              se_naive = sd / sqrt(n_folds), .groups = "drop")
reference <- long |> filter(model == "Clinical-only Cox") |>
  select(fold, metric, reference = value)
differences <- long |> filter(model != "Clinical-only Cox") |>
  left_join(reference, by = c("fold", "metric")) |>
  mutate(improvement = if_else(metric == "uno_cindex", value - reference, reference - value))
delta <- differences |> group_by(model, metric) |>
  summarise(mean = mean(improvement), sd = sd(improvement), n_folds = n(),
              se_naive = sd / sqrt(n_folds), .groups = "drop")
prepare_plot <- function(data) {
  data |> mutate(model = factor(model, levels = model_order),
                   metric = factor(metric, levels = metric_names, labels = metric_labels),
                   kind = case_when(model == "Clinical-only Cox" ~ "Clinical",
                                     model == "Stacked ridge" ~ "Stacked", TRUE ~ "Pathway"))
}
plot_theme <- theme_bw(base_size = 15) +
  theme(axis.text.y = element_text(size = 13), axis.text.x = element_text(size = 12), strip.text = element_text(size = 12),
          plot.margin = margin(10, 14, 10, 10), legend.position = "bottom")
colours <- scale_colour_manual(values = c(Clinical = "firebrick", Stacked = "#1769aa", Pathway = "black"),
                                  guide = "none")
performance_plot <- prepare_plot(summary_metrics) |>
  ggplot(aes(mean, model, colour = kind)) +
  geom_segment(aes(x = mean - se_naive, xend = mean + se_naive, yend = model), linewidth = 0.6) +
  geom_point(size = 2.2) + colours +
  facet_wrap(~ metric, scales = "free_x", nrow = 1) +
  scale_y_discrete(labels = model_label) + scale_x_continuous(n.breaks = 3) +
  labs(title = "Five-fold held-out performance", x = "Metric value", y = NULL) + plot_theme
curve_summary <- a$curves |> group_by(model, time) |>
  summarise(mean = mean(brier_score), se_naive = sd(brier_score) / sqrt(n()), .groups = "drop")
brier_plot <- curve_summary |>
  ggplot(aes(time, mean, colour = model, fill = model, linetype = model)) +
  geom_ribbon(aes(ymin = pmax(0, mean - se_naive), ymax = mean + se_naive),
                alpha = 0.13, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = c("Clinical-only Cox" = "firebrick", "Stacked ridge" = "#1769aa")) +
  scale_fill_manual(values = c("Clinical-only Cox" = "firebrick", "Stacked ridge" = "#1769aa"), guide = "none") +
  labs(title = "Held-out survival-probability error", x = "Follow-up time (days)",
         y = "Mean IPCW Brier score", colour = "Model", linetype = "Model") + plot_theme

# Absolute correlations avoid cancellation when PC1 signs reverse across folds.
correlations <- map_dfr(seq_len(a$settings$k), function(fold) {
  f <- readRDS(file.path(input_folder, paste0("nested_fold_", fold, ".rds")))
  imap(f$models, function(result, pathway) {
    x <- result$training_scores
    tibble(fold, model = pathway,
             comparison = c("PC1 / GSVA", "PC1 / ssGSEA", "GSVA / ssGSEA"),
             abs_r = abs(c(cor(x$PC1, x$GSVA), cor(x$PC1, x$ssGSEA), cor(x$GSVA, x$ssGSEA))))
  }) |> list_rbind()
})
correlation_summary <- correlations |> group_by(model, comparison) |>
  summarise(mean_abs_r = mean(abs_r), .groups = "drop") |>
  mutate(model = factor(model, levels = fixed_pathways))
correlation_plot <- correlation_summary |>
  ggplot(aes(comparison, model, fill = mean_abs_r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", mean_abs_r)), size = 5) +
  scale_fill_gradient(low = "white", high = "#1769aa", limits = c(0, 1),
                        breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"),
                        guide = guide_colourbar(barwidth = grid::unit(2.5, "in"),
                                               barheight = grid::unit(0.15, "in"),
                                               title.position = "top")) +
  scale_y_discrete(labels = model_label) +
  labs(title = "Overlap between pathway scores", x = NULL, y = NULL, fill = "Mean |r|") +
  plot_theme + theme(panel.grid = element_blank())

plots <- list("01_cv_performance" = performance_plot,
                "03_cv_brier_curves" = brier_plot, "04_score_correlations" = correlation_plot)
for (name in names(plots)) {
  for (extension in c("pdf", "png")) {
    ggsave(file.path(figure_folder, paste0(name, ".", extension)), plots[[name]],
             width = 8.5, height = if (name == "03_cv_brier_curves") 6 else 9, dpi = 300)
  }
}

# Table 1: cohort flow, clinical characteristics, and reverse-KM follow-up.
cohort <- a$folds
excluded <- read_csv(file.path(input_folder, "excluded_patients.csv"), show_col_types = FALSE)
median_iqr <- function(x) sprintf("%.1f (%.1f–%.1f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
count_percent <- function(n) sprintf("%d (%.1f%%)", n, 100 * n / nrow(cohort))
reverse_km <- survfit(Surv(OS_time, 1 - OS_event) ~ 1, data = cohort)
followup <- unname(summary(reverse_km)$table["median"])
cohort_table <- tibble(
  Characteristic = c("Patients with valid survival data", "Excluded: missing/invalid clinical covariates",
                       "Patients included", "Observed deaths, n (%)", "Censored, n (%)",
                       "Age, median (IQR), years", "Median follow-up (reverse KM), days",
                       "Female, n (%)", "Male, n (%)", "Stage I, n (%)", "Stage II, n (%)",
                       "Stage III, n (%)", "Stage IV, n (%)"),
  Value = c(nrow(cohort) + nrow(excluded), nrow(excluded), nrow(cohort),
              count_percent(sum(cohort$OS_event)), count_percent(sum(cohort$OS_event == 0)),
              median_iqr(cohort$age), if (is.finite(followup)) sprintf("%.1f", followup) else "Not reached",
              count_percent(sum(cohort$sex == "female")), count_percent(sum(cohort$sex == "male")),
              map_chr(c("stagei", "stageii", "stageiii", "stageiv"),
                        ~ count_percent(sum(cohort$stage == .x))))
)
# Table 2: compact model IDs correspond to the full names in the plots.
performance_table <- summary_metrics |>
  mutate(Model = if_else(model %in% fixed_pathways, unname(ids[model]), model),
           estimate = sprintf("%.3f ± %.3f", mean, se_naive)) |>
  select(Model, metric, estimate) |>
  pivot_wider(names_from = metric, values_from = estimate) |>
  arrange(match(Model, c("Clinical-only Cox", "Stacked ridge", unname(ids)))) |>
  rename(`C-index` = uno_cindex, Brier = brier_score, IBS = integrated_brier_score) |>
  select(Model, `C-index`, Brier, IBS)

write_csv(metrics, file.path(figure_folder, "fold_metrics.csv"))
write_csv(summary_metrics, file.path(figure_folder, "metric_summary.csv"))
write_csv(delta, file.path(figure_folder, "paired_improvement_summary.csv"))
write_csv(curve_summary, file.path(report_folder, "brier_curve_summary.csv"))
write_csv(correlations, file.path(report_folder, "training_fold_correlations.csv"))
write_csv(correlation_summary, file.path(report_folder, "score_correlation_summary.csv"))
write_csv(cohort_table, file.path(report_folder, "patient_characteristics.csv"))
write_csv(performance_table, file.path(report_folder, "performance_table.csv"))
write_csv(tibble(id = ids, pathway = fixed_pathways), file.path(report_folder, "pathway_key.csv"))
write_csv(tibble(outer_fold = seq_len(a$settings$k), penalty = a$selected_penalties),
            file.path(report_folder, "selected_penalties.csv"))

# Draft prose is separate from thesis.qmd: rerunning never overwrites your edits.
number <- function(model, metric) summary_metrics |> filter(.data$model == .env$model, .data$metric == .env$metric) |> pull(mean)
se <- function(model, metric) summary_metrics |> filter(.data$model == .env$model, .data$metric == .env$metric) |> pull(se_naive)
estimate <- function(model, metric) sprintf("%.3f (SE %.3f)", number(model, metric), se(model, metric))
clinical <- map_chr(metric_names, ~ estimate("Clinical-only Cox", .x))
ridge <- map_chr(metric_names, ~ estimate("Stacked ridge", .x))
ridge_delta <- delta |> filter(model == "Stacked ridge")
metric_range <- summary_metrics |> filter(model %in% fixed_pathways, metric == "uno_cindex") |> pull(mean) |> range()
mean_pathway_cor <- correlation_summary |> filter(comparison == "GSVA / ssGSEA") |> pull(mean_abs_r) |> range()
penalties <- paste(a$selected_penalties, collapse = ", ")
changes <- ridge_delta$mean[match(metric_names, ridge_delta$metric)]
assessment <- if (all(changes > 0)) {
  "Stacking had more favourable numerical means for all three metrics."
} else if (all(changes < 0)) {
  "Clinical-only Cox had more favourable numerical means for all three metrics."
} else {
  "The direction of the differences varied across metrics; stacking did not show uniform numerical improvement."
}
curve_wide <- curve_summary |> select(model, time, mean) |>
  pivot_wider(names_from = model, values_from = mean)
fraction_higher <- mean(curve_wide[["Stacked ridge"]] > curve_wide[["Clinical-only Cox"]])
curve_assessment <- if (fraction_higher > 0.5) {
  "The mean stacked-ridge curve was higher (worse) than clinical-only Cox over most of the evaluated grid."
} else {
  "The mean stacked-ridge curve was lower than or equal to clinical-only Cox over at least half of the evaluated grid."
}
md_table <- function(x) paste(knitr::kable(x, format = "pipe", align = c("l", rep("r", ncol(x)-1))), collapse = "\n")
body <- c(
  "## Results", "", "### Study cohort", "",
  sprintf("Of %d patients with valid survival outcomes, %d were excluded for missing or invalid clinical covariates, leaving %d patients, including %d observed deaths. The clinical characteristics are summarized in @tbl-cohort. Each included patient was evaluated once in an outer test fold; the same fold assignments were used for every predictive model.",
            nrow(cohort)+nrow(excluded), nrow(excluded), nrow(cohort), sum(cohort$OS_event)), "",
  md_table(cohort_table), "", ": Patient characteristics for the complete-case analysis. Percentages use the included cohort as denominator; follow-up is estimated by reverse Kaplan–Meier. {#tbl-cohort}", "",
  "```{=latex}", "\\newpage", "```", "",
  "### Cross-validated predictive performance", "",
  sprintf("The clinical-only Cox model had a mean C-index of %s, Brier score of %s, and IBS of %s. The corresponding stacked-ridge estimates were %s, %s, and %s (@tbl-cv-performance; @fig-pathway-performance). The ten individual pathway models had mean C-indices ranging from %.3f to %.3f. These comparisons are descriptive and do not establish statistical superiority.",
            clinical[1], clinical[2], clinical[3], ridge[1], ridge[2], ridge[3], metric_range[1], metric_range[2]), "",
  md_table(performance_table), "", ": Five-fold held-out metrics, reported as mean ± naive fold SE (SD divided by sqrt(5)). P01–P10 identify pathways labelled in the figures. Brier score and C-index are evaluated at 809 days; IBS covers 219.2–809 days. {#tbl-cv-performance}", "",
  "![Held-out performance across five common outer folds. Clinical-only Cox is red, stacked ridge is blue, and individual pathway models are black. Points show fold means; bars show ± one naive fold-based SE. Overlapping training sets make these descriptive SEs approximate; they are not confidence intervals.](figures/reactome_cv/01_cv_performance.pdf){#fig-pathway-performance width=\"100%\" fig-align=\"center\" fig-pos=\"H\"}", "",
  "### Added value beyond clinical information", "",
  sprintf("The stacked model's mean paired changes relative to clinical-only Cox were %.3f for C-index improvement, %.3f for Brier-score improvement, and %.3f for IBS improvement. Positive values favour stacking; negative values favour the clinical reference. Differences were calculated within each test fold before averaging. These descriptive differences are not formal significance tests. The signs and magnitudes should be considered jointly; a gain in discrimination alone does not establish better survival-probability prediction.",
            ridge_delta$mean[match(metric_names, ridge_delta$metric)][1],
            ridge_delta$mean[match(metric_names, ridge_delta$metric)][2],
            ridge_delta$mean[match(metric_names, ridge_delta$metric)][3]), "", assessment, "",
  "```{=latex}", "\\Needspace{0.70\\textheight}", "```", "",
  "### Prediction error over follow-up", "",
  "The cross-validated Brier curves compare the clinical-only and stacked models throughout the common evaluation interval (@fig-cv-brier). A lower curve indicates lower probability-prediction error at that time. The IBS summarizes this error over the full interval; a difference at one time should not be interpreted as evidence of improvement throughout follow-up.", "", curve_assessment, "",
  "![Mean held-out IPCW Brier curves across the five outer folds, with shaded bands showing ± one naive fold SE. Both models use the same test patients, time grid and censoring-weight convention. Lower values are better. Bands are descriptive and are not confidence intervals or significance tests.](figures/reactome_cv/03_cv_brier_curves.pdf){#fig-cv-brier width=\"100%\" fig-align=\"center\" fig-pos=\"H\"}", "",
  "```{=latex}", "\\Needspace{0.90\\textheight}", "```", "",
  "### Overlap between pathway scores", "",
  sprintf("The mean absolute correlations between GSVA and ssGSEA ranged from %.2f to %.2f across the ten pathways (@fig-score-correlations), indicating substantial overlap in their expression information. Absolute correlations are used because the sign of PC1 can reverse between fitted folds. This plot describes score redundancy, not survival association or prediction accuracy.", mean_pathway_cor[1], mean_pathway_cor[2]), "",
  "![Mean absolute Pearson correlations across the five outer training partitions for each pathway. Values near one indicate stronger linear redundancy. Absolute correlations avoid cancellation from arbitrary PC1 sign reversals. Overlapping training sets are not independent replications, and no survival outcomes enter the correlations.](figures/reactome_cv/04_score_correlations.pdf){#fig-score-correlations width=\"100%\" fig-align=\"center\" fig-pos=\"H\"}", "",
  "### Interpretation and scope", "",
  sprintf("The stacked ridge penalty was selected independently inside each outer training set; the five selected penalties were %s. This pilot compares clinical-only, ten pathway-specific, and stacked prediction on shared folds. It does not yet establish the relative merits of gene-level models, single versus vector pathway scores, or biological versus matched random gene sets. The pathway identities and prediction window were fixed in the earlier exploratory analysis, so the present results are internal exploratory validation rather than independent confirmation. GSVA remains cohort-dependent because held-out batches are scored together.", penalties), ""
)
write_lines(body, file.path(report_folder, "results_draft.qmd"))
message("Three figures, two tables, and editable draft saved. Copy the draft into thesis.qmd only when ready.")
