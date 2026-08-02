#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  normalizePath("analysis")
}
repo_dir <- normalizePath(file.path(script_dir, ".."))
source_dir <- file.path(repo_dir, "results", "intermediate", "full_231_significance")
binary_dir <- file.path(repo_dir, "results", "binary_screening")
broad_dir <- file.path(repo_dir, "results", "intermediate", "broad_vs_specific")
output_dir <- file.path(repo_dir, "results", "publication")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

app_levels <- c("AI Skin Scanner", "Model Dermatol", "ChatGPT", "Claude")
metric_levels <- c("Top-1 diagnosis", "Top-3 diagnosis", "Benign/malignant triage")
tone_levels <- c("I-II", "III-IV", "V-VI")
disease_levels <- c(
  "basal cell carcinoma",
  "benign nevus",
  "melanoma spectrum",
  "squamous cell carcinoma"
)

app_labels <- c(
  "AI Skin Scanner" = "AI Skin\nScanner",
  "Model Dermatol" = "Model\nDermatol",
  "ChatGPT" = "ChatGPT",
  "Claude" = "Claude"
)

disease_labels <- c(
  "basal cell carcinoma" = "Basal cell carcinoma",
  "benign nevus" = "Benign nevus",
  "melanoma spectrum" = "Melanoma spectrum",
  "squamous cell carcinoma" = "Squamous cell carcinoma"
)

app_palette <- c(
  "AI Skin Scanner" = "#D55E00",
  "Model Dermatol" = "#0072B2",
  "ChatGPT" = "#009E73",
  "Claude" = "#CC79A7"
)

family_palette <- c(
  "Smartphone apps" = "#B85C38",
  "LLMs" = "#2A9D8F"
)

base_theme <- theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(
    plot.title = element_text(size = 12, face = "bold", lineheight = 1.0),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    panel.spacing = unit(0.7, "lines"),
    plot.margin = margin(5, 8, 5, 5)
  )

save_figure <- function(plot, filename, width, height) {
  ggsave(
    filename = file.path(output_dir, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  ggsave(
    filename = file.path(output_dir, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "pdf",
    bg = "white"
  )
}

overall <- read_csv(
  file.path(source_dir, "overall_by_app_with_significance_full231.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    metric = factor(metric, levels = metric_levels),
    app = factor(app, levels = rev(app_levels)),
    value_label = percent(accuracy, accuracy = 0.1)
  )

omnibus <- read_csv(
  file.path(source_dir, "overall_omnibus_tests_full231.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    metric = factor(metric, levels = metric_levels),
    label = if_else(
      cochran_q_p < 0.001,
      "Overall P<.001",
      paste0("Overall P=", sub("^0", "", sprintf("%.3f", cochran_q_p)))
    )
  )

fig1 <- ggplot(overall, aes(x = accuracy, y = app, color = app)) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0.16,
    linewidth = 0.45
  ) +
  geom_point(size = 2.4) +
  geom_text(
    aes(label = value_label),
    hjust = -0.2,
    size = 2.5,
    color = "black",
    family = "Helvetica"
  ) +
  geom_text(
    data = omnibus,
    aes(x = 0.97, y = 4.55, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    size = 2.5,
    family = "Helvetica"
  ) +
  facet_wrap(~ metric, nrow = 1) +
  scale_color_manual(values = app_palette, guide = "none") +
  scale_y_discrete(labels = app_labels) +
  scale_x_continuous(
    limits = c(0.2, 1),
    breaks = seq(0.25, 1, 0.25),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.09))
  ) +
  labs(x = "Accuracy (95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  base_theme +
  theme(
    panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

save_figure(fig1, "figure1_overall_application_performance", 7.2, 3.25)

family_values <- read_csv(
  file.path(source_dir, "family_metric_plot_values_full231.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    metric = factor(metric, levels = metric_levels),
    family = factor(family, levels = c("Smartphone apps", "LLMs")),
    value_label = percent(mean, accuracy = 0.1)
  )

family_tests <- read_csv(
  file.path(source_dir, "family_metric_significance_full231.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    metric = factor(metric, levels = metric_levels),
    p_label = paste0("Holm P=", sub("^0", "", sprintf("%.3f", p_holm)))
  )

fig2 <- ggplot(family_values, aes(x = mean, y = family, color = family)) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0.15,
    linewidth = 0.5
  ) +
  geom_point(size = 2.6) +
  geom_text(
    aes(label = value_label),
    hjust = -0.2,
    size = 2.6,
    color = "black",
    family = "Helvetica"
  ) +
  geom_text(
    data = family_tests,
    aes(x = 0.97, y = 2.35, label = p_label),
    inherit.aes = FALSE,
    hjust = 1,
    size = 2.5,
    family = "Helvetica"
  ) +
  facet_wrap(~ metric, nrow = 1) +
  scale_color_manual(values = family_palette, guide = "none") +
  scale_x_continuous(
    limits = c(0.3, 1),
    breaks = seq(0.4, 1, 0.2),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.09))
  ) +
  labs(x = "Mean accuracy across the two applications in each family (95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  base_theme +
  theme(
    panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

save_figure(fig2, "figure2_smartphone_vs_llm_performance", 7.2, 2.85)

disease <- read_csv(
  file.path(source_dir, "disease_top1_by_app_full231.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    disease_group = factor(disease_group, levels = disease_levels),
    app = factor(app, levels = rev(app_levels)),
    value_label = percent(accuracy, accuracy = 0.1)
  )

fig3 <- ggplot(disease, aes(x = accuracy, y = app, color = app)) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0.16,
    linewidth = 0.45
  ) +
  geom_point(size = 2.4) +
  geom_text(
    aes(label = value_label),
    hjust = -0.18,
    size = 2.45,
    color = "black",
    family = "Helvetica"
  ) +
  facet_wrap(
    ~ disease_group,
    ncol = 2,
    labeller = as_labeller(disease_labels)
  ) +
  scale_color_manual(values = app_palette, guide = "none") +
  scale_y_discrete(labels = app_labels) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.09))
  ) +
  labs(x = "Top-1 diagnostic accuracy (95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  base_theme +
  theme(
    panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

save_figure(fig3, "figure3_disease_specific_top1_performance", 7.2, 5.0)

tone <- read_csv(
  file.path(source_dir, "tone_top1_by_app_full231.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    app = factor(app, levels = app_levels),
    tone_group = factor(tone_group, levels = tone_levels)
  )

tone_annotation <- tibble(
  app = factor("Claude", levels = app_levels),
  tone_group = factor("III-IV", levels = tone_levels),
  x = 2,
  y = 0.89,
  label = "Overall Holm P=.039\nOrdinal Holm P=.012"
)

fig4 <- ggplot(tone, aes(x = tone_group, y = accuracy, group = 1, color = app)) +
  geom_line(linewidth = 0.55) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12, linewidth = 0.45) +
  geom_point(size = 2.3) +
  geom_text(
    data = tone_annotation,
    aes(x = tone_group, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.5,
    family = "Helvetica"
  ) +
  facet_wrap(~ app, nrow = 1, labeller = as_labeller(app_labels)) +
  scale_color_manual(values = app_palette, guide = "none") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(x = "Grouped Fitzpatrick skin type", y = "Top-1 diagnostic accuracy (95% CI)") +
  base_theme +
  theme(
    panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.35),
    axis.text.x = element_text(size = 7.5)
  )

save_figure(fig4, "figure4_skin_tone_stratified_top1_performance", 7.2, 3.45)

confusion <- read_csv(
  file.path(binary_dir, "binary_confusion_metrics_overall_and_by_tone.csv"),
  show_col_types = FALSE
) %>%
  filter(target == "malignant_vs_benign", stratum == "Overall")

confusion_long <- bind_rows(
  confusion %>% transmute(app, metric = "Sensitivity", estimate = sensitivity, low = sensitivity_ci_low, high = sensitivity_ci_high),
  confusion %>% transmute(app, metric = "Specificity", estimate = specificity, low = specificity_ci_low, high = specificity_ci_high),
  confusion %>% transmute(app, metric = "PPV", estimate = ppv, low = ppv_ci_low, high = ppv_ci_high),
  confusion %>% transmute(app, metric = "NPV", estimate = npv, low = npv_ci_low, high = npv_ci_high)
)

auc_metrics <- function(truth, score) {
  keep <- !is.na(truth) & !is.na(score)
  truth <- as.integer(truth[keep])
  score <- score[keep]
  n_pos <- sum(truth == 1)
  n_neg <- sum(truth == 0)
  ranks <- rank(score, ties.method = "average")
  auroc <- (sum(ranks[truth == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

  thresholds <- sort(unique(score), decreasing = TRUE)
  pr <- map_dfr(thresholds, function(thr) {
    pred <- score >= thr
    tp <- sum(pred & truth == 1)
    fp <- sum(pred & truth == 0)
    tibble(
      recall = tp / n_pos,
      precision = if_else(tp + fp > 0, tp / (tp + fp), 1)
    )
  }) %>%
    arrange(recall, desc(precision)) %>%
    distinct(recall, .keep_all = TRUE) %>%
    bind_rows(tibble(recall = 0, precision = 1), .) %>%
    arrange(recall)

  auprc <- sum(pr$precision[-1] * diff(pr$recall))
  c(AUROC = auroc, AUPRC = auprc)
}

case_level <- read_csv(
  file.path(binary_dir, "case_level_binary_predictions_and_scores.csv"),
  show_col_types = FALSE
) %>%
  filter(target == "malignant_vs_benign")

set.seed(20260725)
bootstrap_auc <- function(df, replicates = 2000) {
  positive <- which(df$truth_positive)
  negative <- which(!df$truth_positive)
  observed <- auc_metrics(df$truth_positive, df$score)

  boots <- replicate(replicates, {
    sampled <- c(
      sample(positive, length(positive), replace = TRUE),
      sample(negative, length(negative), replace = TRUE)
    )
    auc_metrics(df$truth_positive[sampled], df$score[sampled])
  })

  tibble(
    metric = names(observed),
    estimate = as.numeric(observed),
    low = apply(boots, 1, quantile, probs = 0.025, na.rm = TRUE),
    high = apply(boots, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

discrimination_ci_path <- file.path(
  output_dir,
  "malignant_vs_benign_discrimination_bootstrap_ci.csv"
)

discrimination_ci <- case_level %>%
  group_by(app) %>%
  group_modify(~ bootstrap_auc(.x)) %>%
  ungroup()

write_csv(discrimination_ci, discrimination_ci_path)

screening <- bind_rows(confusion_long, discrimination_ci) %>%
  mutate(
    metric = factor(metric, levels = c("Sensitivity", "Specificity", "PPV", "NPV", "AUROC", "AUPRC")),
    app = factor(app, levels = rev(app_levels)),
    value_label = sprintf("%.2f", estimate)
  )

fig5 <- ggplot(screening, aes(x = estimate, y = app, color = app)) +
  geom_errorbar(
    aes(xmin = low, xmax = high),
    orientation = "y",
    width = 0.16,
    linewidth = 0.45
  ) +
  geom_point(size = 2.35) +
  geom_text(
    aes(label = value_label),
    hjust = -0.25,
    size = 2.35,
    color = "black",
    family = "Helvetica"
  ) +
  facet_wrap(~ metric, ncol = 3) +
  scale_color_manual(values = app_palette, guide = "none") +
  scale_y_discrete(labels = app_labels) +
  scale_x_continuous(
    limits = c(0.3, 1),
    breaks = seq(0.4, 1, 0.2),
    labels = number_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.1))
  ) +
  labs(x = "Estimate (95% CI)", y = NULL) +
  coord_cartesian(clip = "off") +
  base_theme +
  theme(
    panel.grid.major.x = element_line(color = "#E6E6E6", linewidth = 0.35),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )

save_figure(fig5, "figure5_malignancy_screening_performance", 7.2, 5.25)

broad <- read_csv(
  file.path(broad_dir, "broad_vs_specific_top1_by_disease_app.csv"),
  show_col_types = FALSE
) %>%
  group_by(app, mapping) %>%
  summarise(
    n = sum(n),
    correct_n = sum(correct_n),
    accuracy = correct_n / n,
    .groups = "drop"
  ) %>%
  mutate(
    app = factor(app, levels = app_levels),
    mapping = factor(mapping, levels = c("Broad allowed", "Specific only"))
  )

fig_s1 <- ggplot(broad, aes(x = app, y = accuracy, fill = mapping)) +
  geom_col(
    position = position_dodge(width = 0.74),
    width = 0.64
  ) +
  geom_text(
    aes(label = percent(accuracy, accuracy = 0.1)),
    position = position_dodge(width = 0.74),
    vjust = -0.45,
    size = 2.9,
    color = "black",
    family = "Helvetica"
  ) +
  scale_fill_manual(
    values = c("Broad allowed" = "#4C78A8", "Specific only" = "#F58518")
  ) +
  scale_x_discrete(labels = app_labels) +
  scale_y_continuous(
    limits = c(0, 0.72),
    breaks = seq(0, 0.6, 0.2),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Broad terminology changes only Model Dermatol's apparent top-1 accuracy",
    subtitle = "Broad mapping accepts category-level terms; specific mapping requires a disease-level diagnosis",
    x = NULL,
    y = "Top-1 diagnostic accuracy",
    fill = NULL
  ) +
  base_theme +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(color = "#E6E6E6", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 8)
  )

save_figure(fig_s1, "supplementary_figure1_diagnosis_specificity_sensitivity", 7.2, 4.6)

write_csv(overall, file.path(output_dir, "overall_application_performance.csv"))
write_csv(family_values, file.path(output_dir, "application_family_performance.csv"))
write_csv(disease, file.path(output_dir, "disease_specific_top1_performance.csv"))
write_csv(tone, file.path(output_dir, "skin_tone_stratified_top1_performance.csv"))
write_csv(screening, file.path(output_dir, "malignancy_screening_performance.csv"))

message("Publication results figures written to: ", output_dir)
