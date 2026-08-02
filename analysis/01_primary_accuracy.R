#!/usr/bin/env Rscript

# Analysis script for Trinidad Lab fitzpatrick dataset paper figures/stats.
# It reads the case-level CSV, reshapes app results to long format, writes
# summary tables, runs core statistical tests, and saves publication-style plots.

# If needed, install packages first with:
# install.packages(c("dplyr", "tidyr", "readr", "stringr", "forcats",
#                    "ggplot2", "purrr", "broom", "scales", "patchwork"))

required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "forcats",
  "ggplot2", "purrr", "broom", "scales", "patchwork"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them with install.packages(...) and re-run the script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(purrr)
  library(broom)
  library(scales)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  normalizePath("analysis")
}
repo_dir <- normalizePath(file.path(script_dir, ".."))
input_path <- Sys.getenv(
  "INPUT_PATH",
  unset = file.path(repo_dir, "data", "fitzpatrick", "fitzpatrick_analysis_231.csv")
)
output_dir <- Sys.getenv(
  "OUTPUT_DIR",
  unset = file.path(repo_dir, "results", "primary_accuracy")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
disease_output_dir <- file.path(output_dir, "disease_specific")
dir.create(disease_output_dir, recursive = TRUE, showWarnings = FALSE)

apps <- c(
  a = "AI Skin Scanner",
  m = "Model Dermatol",
  c = "ChatGPT",
  cl = "Claude"
)

outcomes <- c(
  dx_1_correct = "Top-1 diagnosis accuracy",
  any_dx_correct = "Any diagnosis accuracy",
  dx_1_malignancy_correct = "Benign/malignant accuracy"
)

fitz_group_3 <- function(x) {
  case_when(
    x %in% c("1", "2") ~ "I-II",
    x %in% c("3", "4") ~ "III-IV",
    x %in% c("5", "6") ~ "V-VI",
    TRUE ~ "Unknown"
  )
}

slugify <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

wilson_ci <- function(successes, total, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  valid <- !is.na(total) & total > 0
  phat <- rep(NA_real_, length(total))
  denom <- rep(NA_real_, length(total))
  center <- rep(NA_real_, length(total))
  half_width <- rep(NA_real_, length(total))

  phat[valid] <- successes[valid] / total[valid]
  denom[valid] <- 1 + (z^2 / total[valid])
  center[valid] <- (phat[valid] + (z^2 / (2 * total[valid]))) / denom[valid]
  half_width[valid] <- (z / denom[valid]) * sqrt(
    (phat[valid] * (1 - phat[valid]) / total[valid]) + (z^2 / (4 * total[valid]^2))
  )

  ci_low <- rep(NA_real_, length(total))
  ci_high <- rep(NA_real_, length(total))
  ci_low[valid] <- pmax(0, center[valid] - half_width[valid])
  ci_high[valid] <- pmin(1, center[valid] + half_width[valid])

  tibble(
    ci_low = ci_low,
    ci_high = ci_high
  )
}

safe_test_name <- function(tab) {
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (any(expected < 5)) "fisher" else "chisq"
}

run_skin_tone_test <- function(df, outcome_name) {
  tab <- xtabs(correct ~ fitzpatrick_scale + app, data = df %>% filter(outcome == outcome_name))
  # Not used directly; retained as a helper if you want a global test later.
  tab
}

message("Reading data from: ", input_path)

raw <- read_csv(input_path, show_col_types = FALSE) %>%
  mutate(
    qc_flag = if_else(is.na(qc) | qc == "", "No QC flag", qc),
    fitzpatrick_scale = as.character(fitzpatrick_scale),
    fitzpatrick_centaur = as.character(fitzpatrick_centaur),
    fitzpatrick_group_3 = fitz_group_3(fitzpatrick_scale),
    label = str_squish(label),
    label = case_when(
      label %in% c("congenital nevus", "halo nevus") ~ "benign nevus",
      label %in% c("malignant melanoma", "melanoma", "lentigo maligna") ~ "melanoma spectrum",
      TRUE ~ label
    )
  )

long <- raw %>%
  pivot_longer(
    cols = matches("^(a|m|c|cl)_(dx_1_correct|any_dx_correct|dx_1_malignancy_correct)$"),
    names_to = c("app_code", "outcome"),
    names_pattern = "^(a|m|c|cl)_(dx_1_correct|any_dx_correct|dx_1_malignancy_correct)$",
    values_to = "correct"
  ) %>%
  mutate(
    app = factor(unname(apps[app_code]), levels = unname(apps)),
    outcome = factor(unname(outcomes[outcome]), levels = unname(outcomes)),
    correct = as.integer(correct),
    fitzpatrick_centaur = factor(fitzpatrick_centaur, levels = as.character(1:6)),
    fitzpatrick_scale = factor(fitzpatrick_scale, levels = as.character(1:6)),
    fitzpatrick_group_3 = factor(fitzpatrick_group_3, levels = c("I-II", "III-IV", "V-VI")),
    label = fct_infreq(label)
  )

write_csv(long, file.path(output_dir, "fitzpatrick_long_format.csv"))

# Summary tables -------------------------------------------------------------

overall_summary <- long %>%
  group_by(app, outcome) %>%
  summarise(
    n = n(),
    correct_n = sum(correct),
    accuracy = correct_n / n,
    .groups = "drop"
  ) %>%
  bind_cols(wilson_ci(.$correct_n, .$n))

write_csv(overall_summary, file.path(output_dir, "overall_accuracy_by_app.csv"))

skin_tone_summary <- long %>%
  group_by(app, outcome, fitzpatrick_scale) %>%
  summarise(
    n = n(),
    correct_n = sum(correct),
    accuracy = correct_n / n,
    .groups = "drop"
  ) %>%
  bind_cols(wilson_ci(.$correct_n, .$n))

write_csv(skin_tone_summary, file.path(output_dir, "accuracy_by_app_and_scale_tone.csv"))

skin_group_summary <- long %>%
  group_by(app, outcome, fitzpatrick_group_3) %>%
  summarise(
    n = n(),
    correct_n = sum(correct),
    accuracy = correct_n / n,
    .groups = "drop"
  ) %>%
  bind_cols(wilson_ci(.$correct_n, .$n))

write_csv(skin_group_summary, file.path(output_dir, "accuracy_by_app_and_3_tone_group.csv"))

disease_summary <- long %>%
  group_by(app, outcome, label) %>%
  summarise(
    n = n(),
    correct_n = sum(correct),
    accuracy = correct_n / n,
    .groups = "drop"
  ) %>%
  bind_cols(wilson_ci(.$correct_n, .$n))

write_csv(disease_summary, file.path(output_dir, "accuracy_by_app_and_disease.csv"))

disease_tone_summary <- long %>%
  group_by(app, outcome, label, fitzpatrick_scale) %>%
  summarise(
    n = n(),
    correct_n = sum(correct),
    accuracy = correct_n / n,
    .groups = "drop"
  ) %>%
  bind_cols(wilson_ci(.$correct_n, .$n))

write_csv(disease_tone_summary, file.path(output_dir, "accuracy_by_app_disease_and_scale_tone.csv"))

# Statistical tests ----------------------------------------------------------

# 1. Across skin tones within each app: chi-square or Fisher exact on correct vs incorrect.
skin_tone_tests <- long %>%
  group_by(app, outcome) %>%
  group_modify(~ {
    tab <- table(.x$fitzpatrick_scale, .x$correct)
    colnames(tab) <- c("incorrect", "correct")
    test_type <- safe_test_name(tab)
    test_obj <- if (test_type == "fisher") fisher.test(tab) else chisq.test(tab)

    tibble(
      test = if_else(test_type == "fisher", "Fisher exact", "Chi-square"),
      statistic = if (is.null(test_obj$statistic)) NA_real_ else unname(test_obj$statistic),
      df = if (is.null(test_obj$parameter)) NA_real_ else unname(test_obj$parameter),
      p_value = test_obj$p.value
    )
  }) %>%
  ungroup() %>%
  arrange(outcome, app)

write_csv(skin_tone_tests, file.path(output_dir, "skin_tone_tests_within_app.csv"))

# 2. Across apps overall: McNemar for paired comparisons.
paired_wide <- raw %>%
  transmute(
    md5hash,
    a_dx_1_correct, m_dx_1_correct, c_dx_1_correct, cl_dx_1_correct,
    a_any_dx_correct, m_any_dx_correct, c_any_dx_correct, cl_any_dx_correct,
    a_dx_1_malignancy_correct, m_dx_1_malignancy_correct, c_dx_1_malignancy_correct, cl_dx_1_malignancy_correct
  )

outcome_codes <- names(outcomes)
app_pairs <- combn(names(apps), 2, simplify = FALSE)

app_pairwise_tests <- map_dfr(outcome_codes, function(metric) {
  map_dfr(app_pairs, function(pair) {
    x <- paired_wide[[paste0(pair[1], "_", metric)]]
    y <- paired_wide[[paste0(pair[2], "_", metric)]]
    tab <- table(factor(x, levels = c(0, 1)), factor(y, levels = c(0, 1)))
    test_obj <- mcnemar.test(tab, correct = TRUE)

    tibble(
      outcome = outcomes[[metric]],
      app_1 = apps[[pair[1]]],
      app_2 = apps[[pair[2]]],
      discordant_10 = tab["1", "0"],
      discordant_01 = tab["0", "1"],
      p_value = test_obj$p.value
    )
  })
})

app_pairwise_tests <- app_pairwise_tests %>%
  group_by(outcome) %>%
  mutate(p_value_holm = p.adjust(p_value, method = "holm")) %>%
  ungroup()

write_csv(app_pairwise_tests, file.path(output_dir, "paired_app_comparisons_mcnemar.csv"))

# 3. Main regression model: adjust for app, skin tone, and disease.
# This treats each app-case prediction as an observation.
glm_results <- map_dfr(outcome_codes, function(metric) {
  df <- long %>%
    filter(as.character(outcome) == outcomes[[metric]]) %>%
    mutate(
      app = relevel(app, ref = "ChatGPT"),
      fitzpatrick_scale = relevel(fitzpatrick_scale, ref = "1")
    )

  fit <- glm(
    correct ~ app + fitzpatrick_scale + label,
    family = binomial(),
    data = df
  )

  tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(outcome = outcomes[[metric]])
})

write_csv(glm_results, file.path(output_dir, "logistic_regression_odds_ratios.csv"))

# 4. Within-disease skin tone tests for diagnosis accuracy only.
disease_tone_tests <- long %>%
  filter(outcome == "Top-1 diagnosis accuracy") %>%
  group_by(app, label) %>%
  group_modify(~ {
    tab <- table(.x$fitzpatrick_scale, .x$correct)
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(tibble(test = NA_character_, statistic = NA_real_, df = NA_real_, p_value = NA_real_))
    }

    test_type <- safe_test_name(tab)
    test_obj <- if (test_type == "fisher") fisher.test(tab) else chisq.test(tab)

    tibble(
      test = if_else(test_type == "fisher", "Fisher exact", "Chi-square"),
      statistic = if (is.null(test_obj$statistic)) NA_real_ else unname(test_obj$statistic),
      df = if (is.null(test_obj$parameter)) NA_real_ else unname(test_obj$parameter),
      p_value = test_obj$p.value
    )
  }) %>%
  ungroup() %>%
  arrange(app, p_value, label)

write_csv(disease_tone_tests, file.path(output_dir, "within_disease_skin_tone_tests_top1_dx.csv"))

# 5. Within-disease app tests across apps for top-1 diagnosis accuracy.
disease_app_tests <- long %>%
  filter(outcome == "Top-1 diagnosis accuracy") %>%
  group_by(label) %>%
  group_modify(~ {
    tab <- table(.x$app, .x$correct)
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(tibble(test = NA_character_, statistic = NA_real_, df = NA_real_, p_value = NA_real_))
    }

    test_type <- safe_test_name(tab)
    test_obj <- if (test_type == "fisher") fisher.test(tab) else chisq.test(tab)

    tibble(
      test = if_else(test_type == "fisher", "Fisher exact", "Chi-square"),
      statistic = if (is.null(test_obj$statistic)) NA_real_ else unname(test_obj$statistic),
      df = if (is.null(test_obj$parameter)) NA_real_ else unname(test_obj$parameter),
      p_value = test_obj$p.value
    )
  }) %>%
  ungroup() %>%
  arrange(p_value, label)

write_csv(disease_app_tests, file.path(output_dir, "within_disease_app_tests_top1_dx.csv"))

# Figures -------------------------------------------------------------------

theme_set(
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8))
    )
)

plot_overall <- ggplot(overall_summary, aes(x = app, y = accuracy, fill = app)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.18) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Overall accuracy by app",
    x = NULL,
    y = "Accuracy"
  )

overall_sig_brackets <- app_pairwise_tests %>%
  filter(!is.na(p_value_holm), p_value_holm < 0.05) %>%
  mutate(
    x1 = match(app_1, levels(overall_summary$app)),
    x2 = match(app_2, levels(overall_summary$app)),
    x_left = pmin(x1, x2),
    x_right = pmax(x1, x2),
    label = case_when(
      p_value_holm < 0.001 ~ "***",
      p_value_holm < 0.01 ~ "**",
      TRUE ~ "*"
    )
  ) %>%
  left_join(
    overall_summary %>%
      group_by(outcome) %>%
      summarise(panel_top = max(ci_high), .groups = "drop"),
    by = "outcome"
  ) %>%
  arrange(outcome, x_right - x_left, x_left, x_right) %>%
  group_by(outcome) %>%
  mutate(
    bracket_row = row_number(),
    y = pmin(0.98, panel_top + 0.05 + (bracket_row - 1) * 0.07),
    y_tip = y - 0.025,
    x_mid = (x_left + x_right) / 2
  ) %>%
  ungroup()

plot_overall <- plot_overall +
  geom_segment(
    data = overall_sig_brackets,
    aes(x = x_left, xend = x_right, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = overall_sig_brackets,
    aes(x = x_left, xend = x_left, y = y_tip, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = overall_sig_brackets,
    aes(x = x_right, xend = x_right, y = y_tip, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = overall_sig_brackets,
    aes(x = x_mid, y = y + 0.015, label = label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 4
  )

ggsave(
  filename = file.path(output_dir, "overall_accuracy_by_app.png"),
  plot = plot_overall,
  width = 12,
  height = 5,
  dpi = 300
)

plot_skin_tone <- ggplot(skin_tone_summary, aes(x = fitzpatrick_scale, y = accuracy, color = app, group = app)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Accuracy by Fitzpatrick skin tone",
    x = "Fitzpatrick skin tone",
    y = "Accuracy",
    color = "App"
  )

ggsave(
  filename = file.path(output_dir, "accuracy_by_skin_tone.png"),
  plot = plot_skin_tone,
  width = 12,
  height = 5,
  dpi = 300
)

heatmap_data <- disease_tone_summary %>%
  filter(outcome == "Top-1 diagnosis accuracy") %>%
  mutate(label = fct_reorder(label, accuracy, .fun = mean, .desc = TRUE))

plot_heatmap <- ggplot(heatmap_data, aes(x = fitzpatrick_scale, y = label, fill = accuracy)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ app, ncol = 2) +
  scale_fill_gradient(low = "#f4efe6", high = "#0f5c5c", labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Top-1 diagnosis accuracy within disease and skin tone",
    x = "Fitzpatrick skin tone",
    y = "Disease",
    fill = "Accuracy"
  ) +
  theme(axis.text.y = element_text(size = 8))

ggsave(
  filename = file.path(output_dir, "heatmap_top1_dx_by_disease_skin_tone.png"),
  plot = plot_heatmap,
  width = 12,
  height = 9,
  dpi = 300
)

plot_grouped <- ggplot(skin_group_summary, aes(x = fitzpatrick_group_3, y = accuracy, fill = app)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    position = position_dodge(width = 0.75),
    width = 0.18
  ) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Accuracy by grouped skin tone",
    x = "Grouped Fitzpatrick tone",
    y = "Accuracy",
    fill = "App"
  )

ggsave(
  filename = file.path(output_dir, "accuracy_by_grouped_skin_tone.png"),
  plot = plot_grouped,
  width = 12,
  height = 5,
  dpi = 300
)

plot_disease_by_app <- disease_summary %>%
  mutate(label_chr = as.character(label)) %>%
  group_by(label_chr) %>%
  mutate(max_accuracy_for_order = max(accuracy, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    label_chr = fct_reorder(label_chr, max_accuracy_for_order, .desc = TRUE)
  ) %>%
  ggplot(aes(x = accuracy, y = label_chr, color = app)) +
  geom_point(size = 2.4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.18) +
  facet_wrap(~ outcome, ncol = 1) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Disease-level accuracy within each app",
    x = "Accuracy",
    y = "Disease",
    color = "App"
  )

ggsave(
  filename = file.path(output_dir, "accuracy_by_disease_within_app.png"),
  plot = plot_disease_by_app,
  width = 12,
  height = 10,
  dpi = 300
)

# Disease-specific outputs ---------------------------------------------------

disease_plot_data <- disease_tone_summary %>%
  filter(outcome %in% c("Top-1 diagnosis accuracy", "Benign/malignant accuracy")) %>%
  mutate(
    label_chr = as.character(label),
    outcome = factor(outcome, levels = c("Top-1 diagnosis accuracy", "Benign/malignant accuracy"))
  )

unique_diseases <- sort(unique(disease_plot_data$label_chr))

walk(unique_diseases, function(disease_name) {
  disease_slug <- slugify(disease_name)

  disease_summary_csv <- disease_plot_data %>%
    filter(label_chr == disease_name) %>%
    select(label = label_chr, app, outcome, fitzpatrick_scale, n, correct_n, accuracy, ci_low, ci_high)

  write_csv(
    disease_summary_csv,
    file.path(disease_output_dir, paste0(disease_slug, "_summary.csv"))
  )

  disease_plot <- ggplot(
    disease_plot_data %>% filter(label_chr == disease_name),
    aes(x = fitzpatrick_scale, y = accuracy, color = app, group = app)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12) +
    facet_wrap(~ outcome, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      title = paste("Accuracy by app and skin tone:", disease_name),
      x = "Fitzpatrick skin tone",
      y = "Accuracy",
      color = "App"
    ) +
    theme(legend.position = "bottom")

  ggsave(
    filename = file.path(disease_output_dir, paste0(disease_slug, "_accuracy_by_tone.png")),
    plot = disease_plot,
    width = 11,
    height = 4.8,
    dpi = 300
  )
})

focus_diseases <- c("squamous cell carcinoma", "basal cell carcinoma", "melanoma", "malignant melanoma")

focus_plot_data <- disease_plot_data %>%
  filter(label_chr %in% focus_diseases)

if (nrow(focus_plot_data) > 0) {
  plot_focus_diseases <- ggplot(
    focus_plot_data,
    aes(x = fitzpatrick_scale, y = accuracy, color = app, group = app)
  ) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
    facet_grid(outcome ~ label_chr) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      title = "Disease-specific accuracy by app and skin tone",
      x = "Fitzpatrick skin tone",
      y = "Accuracy",
      color = "App"
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 0)
    )

  ggsave(
    filename = file.path(output_dir, "focused_disease_accuracy_by_tone.png"),
    plot = plot_focus_diseases,
    width = 14,
    height = 7,
    dpi = 300
  )
}

message("Analysis complete.")
message("Outputs written to: ", output_dir)
