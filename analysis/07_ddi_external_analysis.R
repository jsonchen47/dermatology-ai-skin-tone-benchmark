#!/usr/bin/env Rscript

# DDI analysis workflow for balanced 84-per-tone subset.
# Merges ground truth metadata with ChatGPT and Claude predictions, then writes
# summary tables and publication-style figures for exact diagnosis accuracy and
# benign/malignant accuracy across skin tones and diseases.

required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "forcats",
  "ggplot2", "purrr", "scales"
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
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  normalizePath("analysis")
}
repo_dir <- normalizePath(file.path(script_dir, ".."))
base_dir <- file.path(repo_dir, "data", "ddi")
truth_path <- file.path(base_dir, "ddi_cohort_252.csv")
chat_path <- file.path(base_dir, "chatgpt_outputs_252.csv")
claude_path <- file.path(base_dir, "claude_outputs_252.csv")
output_dir <- file.path(repo_dir, "results", "ddi")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

apps <- c(
  chat = "ChatGPT",
  claude = "Claude"
)

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

  tibble(ci_low = ci_low, ci_high = ci_high)
}

safe_test_name <- function(tab) {
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (any(expected < 5)) "fisher" else "chisq"
}

cochran_armitage_test <- function(successes, totals, scores) {
  total_n <- sum(totals)
  total_success <- sum(successes)
  if (total_n == 0 || total_success == 0 || total_success == total_n) {
    return(tibble(statistic = NA_real_, p_value = NA_real_))
  }

  p_hat <- total_success / total_n
  score_mean <- sum(scores * totals) / total_n
  numerator <- sum(scores * (successes - totals * p_hat))
  variance <- p_hat * (1 - p_hat) * sum(totals * (scores - score_mean)^2)

  if (variance <= 0) {
    return(tibble(statistic = NA_real_, p_value = NA_real_))
  }

  z_value <- numerator / sqrt(variance)

  tibble(
    statistic = z_value,
    p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE)
  )
}

truth <- read_csv(truth_path, show_col_types = FALSE) %>%
  mutate(
    DDI_ID = as.character(DDI_ID),
    true_malignant = str_to_lower(malignant) == "true",
    skin_tone = factor(as.character(skin_tone), levels = c("12", "34", "56"))
  ) %>%
  transmute(
    DDI_ID,
    DDI_file,
    skin_tone,
    disease,
    true_malignant
  )

load_predictions <- function(path, app_code) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(
      DDI_ID = as.character(DDI_ID),
      pred_malignant = str_to_lower(malignant) == "true",
      app = factor(apps[[app_code]], levels = unname(apps))
    ) %>%
    transmute(
      DDI_ID,
      app,
      pred_diagnosis = diagnosis,
      pred_malignant
    )
}

predictions <- bind_rows(
  load_predictions(chat_path, "chat"),
  load_predictions(claude_path, "claude")
)

analysis_long <- truth %>%
  inner_join(predictions, by = "DDI_ID") %>%
  mutate(
    dx_correct = as.integer(pred_diagnosis == disease),
    mal_correct = as.integer(pred_malignant == true_malignant)
  )

write_csv(analysis_long, file.path(output_dir, "ddi_long_analysis.csv"))

overall_summary <- bind_rows(
  analysis_long %>%
    group_by(app) %>%
    summarise(
      outcome = "Top-1 diagnosis accuracy",
      n = n(),
      correct_n = sum(dx_correct),
      accuracy = correct_n / n,
      .groups = "drop"
    ),
  analysis_long %>%
    group_by(app) %>%
    summarise(
      outcome = "Benign/malignant accuracy",
      n = n(),
      correct_n = sum(mal_correct),
      accuracy = correct_n / n,
      .groups = "drop"
    )
) %>%
  bind_cols(wilson_ci(.$correct_n, .$n)) %>%
  mutate(outcome = factor(outcome, levels = c("Top-1 diagnosis accuracy", "Benign/malignant accuracy")))

write_csv(overall_summary, file.path(output_dir, "overall_accuracy_by_app.csv"))

skin_tone_summary <- bind_rows(
  analysis_long %>%
    group_by(app, skin_tone) %>%
    summarise(
      outcome = "Top-1 diagnosis accuracy",
      n = n(),
      correct_n = sum(dx_correct),
      accuracy = correct_n / n,
      .groups = "drop"
    ),
  analysis_long %>%
    group_by(app, skin_tone) %>%
    summarise(
      outcome = "Benign/malignant accuracy",
      n = n(),
      correct_n = sum(mal_correct),
      accuracy = correct_n / n,
      .groups = "drop"
    )
) %>%
  bind_cols(wilson_ci(.$correct_n, .$n)) %>%
  mutate(outcome = factor(outcome, levels = c("Top-1 diagnosis accuracy", "Benign/malignant accuracy")))

write_csv(skin_tone_summary, file.path(output_dir, "accuracy_by_skin_tone.csv"))

disease_summary <- bind_rows(
  analysis_long %>%
    group_by(app, disease) %>%
    summarise(
      outcome = "Top-1 diagnosis accuracy",
      n = n(),
      correct_n = sum(dx_correct),
      accuracy = correct_n / n,
      .groups = "drop"
    ),
  analysis_long %>%
    group_by(app, disease) %>%
    summarise(
      outcome = "Benign/malignant accuracy",
      n = n(),
      correct_n = sum(mal_correct),
      accuracy = correct_n / n,
      .groups = "drop"
    )
) %>%
  bind_cols(wilson_ci(.$correct_n, .$n)) %>%
  mutate(outcome = factor(outcome, levels = c("Top-1 diagnosis accuracy", "Benign/malignant accuracy")))

write_csv(disease_summary, file.path(output_dir, "accuracy_by_disease.csv"))

disease_tone_summary <- bind_rows(
  analysis_long %>%
    group_by(app, disease, skin_tone) %>%
    summarise(
      outcome = "Top-1 diagnosis accuracy",
      n = n(),
      correct_n = sum(dx_correct),
      accuracy = correct_n / n,
      .groups = "drop"
    ),
  analysis_long %>%
    group_by(app, disease, skin_tone) %>%
    summarise(
      outcome = "Benign/malignant accuracy",
      n = n(),
      correct_n = sum(mal_correct),
      accuracy = correct_n / n,
      .groups = "drop"
    )
) %>%
  bind_cols(wilson_ci(.$correct_n, .$n)) %>%
  mutate(outcome = factor(outcome, levels = c("Top-1 diagnosis accuracy", "Benign/malignant accuracy")))

write_csv(disease_tone_summary, file.path(output_dir, "accuracy_by_disease_and_skin_tone.csv"))

paired_app_tests <- bind_rows(
  analysis_long %>%
    select(DDI_ID, app, dx_correct) %>%
    pivot_wider(names_from = app, values_from = dx_correct) %>%
    summarise(
      outcome = "Top-1 diagnosis accuracy",
      discordant_10 = sum(ChatGPT == 1 & Claude == 0),
      discordant_01 = sum(ChatGPT == 0 & Claude == 1),
      p_value = mcnemar.test(table(
        factor(ChatGPT, levels = c(0, 1)),
        factor(Claude, levels = c(0, 1))
      ), correct = TRUE)$p.value
    ),
  analysis_long %>%
    select(DDI_ID, app, mal_correct) %>%
    pivot_wider(names_from = app, values_from = mal_correct) %>%
    summarise(
      outcome = "Benign/malignant accuracy",
      discordant_10 = sum(ChatGPT == 1 & Claude == 0),
      discordant_01 = sum(ChatGPT == 0 & Claude == 1),
      p_value = mcnemar.test(table(
        factor(ChatGPT, levels = c(0, 1)),
        factor(Claude, levels = c(0, 1))
      ), correct = TRUE)$p.value
    )
) %>%
  mutate(
    p_value_holm = p.adjust(p_value, method = "holm"),
    outcome = factor(outcome, levels = c("Top-1 diagnosis accuracy", "Benign/malignant accuracy"))
  )

write_csv(paired_app_tests, file.path(output_dir, "paired_app_comparisons_mcnemar.csv"))

skin_tone_tests <- bind_rows(
  analysis_long %>%
    group_by(app) %>%
    group_modify(~ {
      tab <- table(.x$skin_tone, .x$dx_correct)
      test_type <- safe_test_name(tab)
      test_obj <- if (test_type == "fisher") fisher.test(tab) else chisq.test(tab)
      tibble(
        outcome = "Top-1 diagnosis accuracy",
        test = if_else(test_type == "fisher", "Fisher exact", "Chi-square"),
        statistic = if (is.null(test_obj$statistic)) NA_real_ else unname(test_obj$statistic),
        df = if (is.null(test_obj$parameter)) NA_real_ else unname(test_obj$parameter),
        p_value = test_obj$p.value
      )
    }) %>%
    ungroup(),
  analysis_long %>%
    group_by(app) %>%
    group_modify(~ {
      tab <- table(.x$skin_tone, .x$mal_correct)
      test_type <- safe_test_name(tab)
      test_obj <- if (test_type == "fisher") fisher.test(tab) else chisq.test(tab)
      tibble(
        outcome = "Benign/malignant accuracy",
        test = if_else(test_type == "fisher", "Fisher exact", "Chi-square"),
        statistic = if (is.null(test_obj$statistic)) NA_real_ else unname(test_obj$statistic),
        df = if (is.null(test_obj$parameter)) NA_real_ else unname(test_obj$parameter),
        p_value = test_obj$p.value
      )
    }) %>%
    ungroup()
) %>%
  mutate(outcome = factor(outcome, levels = c("Top-1 diagnosis accuracy", "Benign/malignant accuracy")))

write_csv(skin_tone_tests, file.path(output_dir, "within_app_skin_tone_tests.csv"))

skin_tone_ordinal_tests <- bind_rows(
  analysis_long %>%
    group_by(app) %>%
    group_modify(~ {
      counts <- .x %>%
        group_by(skin_tone) %>%
        summarise(correct_n = sum(dx_correct), n = n(), .groups = "drop") %>%
        arrange(skin_tone)
      trend <- cochran_armitage_test(counts$correct_n, counts$n, c(1, 2, 3))
      tibble(
        outcome = "Top-1 diagnosis accuracy",
        statistic = trend$statistic,
        p_value = trend$p_value
      )
    }) %>%
    ungroup(),
  analysis_long %>%
    group_by(app) %>%
    group_modify(~ {
      counts <- .x %>%
        group_by(skin_tone) %>%
        summarise(correct_n = sum(mal_correct), n = n(), .groups = "drop") %>%
        arrange(skin_tone)
      trend <- cochran_armitage_test(counts$correct_n, counts$n, c(1, 2, 3))
      tibble(
        outcome = "Benign/malignant accuracy",
        statistic = trend$statistic,
        p_value = trend$p_value
      )
    }) %>%
    ungroup()
) %>%
  mutate(outcome = factor(outcome, levels = c("Top-1 diagnosis accuracy", "Benign/malignant accuracy")))

write_csv(skin_tone_ordinal_tests, file.path(output_dir, "within_app_skin_tone_ordinal_trend_tests.csv"))

pair_defs <- list(
  c("12", "34"),
  c("12", "56"),
  c("34", "56")
)

skin_tone_pairwise_tests <- bind_rows(
  analysis_long %>%
    group_by(app) %>%
    group_modify(~ {
      map_dfr(pair_defs, function(pair) {
        sub_df <- .x %>% filter(as.character(skin_tone) %in% pair)
        tab <- table(sub_df$skin_tone, sub_df$dx_correct)
        tibble(
          outcome = "Top-1 diagnosis accuracy",
          comparison = paste(pair, collapse = " vs "),
          p_value = fisher.test(tab)$p.value
        )
      })
    }) %>%
    ungroup(),
  analysis_long %>%
    group_by(app) %>%
    group_modify(~ {
      map_dfr(pair_defs, function(pair) {
        sub_df <- .x %>% filter(as.character(skin_tone) %in% pair)
        tab <- table(sub_df$skin_tone, sub_df$mal_correct)
        tibble(
          outcome = "Benign/malignant accuracy",
          comparison = paste(pair, collapse = " vs "),
          p_value = fisher.test(tab)$p.value
        )
      })
    }) %>%
    ungroup()
) %>%
  group_by(app, outcome) %>%
  mutate(p_value_holm = p.adjust(p_value, method = "holm")) %>%
  ungroup() %>%
  mutate(outcome = factor(outcome, levels = c("Top-1 diagnosis accuracy", "Benign/malignant accuracy")))

write_csv(skin_tone_pairwise_tests, file.path(output_dir, "within_app_skin_tone_pairwise_tests.csv"))

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
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.18) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Overall DDI accuracy by app",
    x = NULL,
    y = "Accuracy"
  )

overall_brackets <- paired_app_tests %>%
  filter(!is.na(p_value_holm), p_value_holm < 0.05) %>%
  left_join(
    overall_summary %>%
      group_by(outcome) %>%
      summarise(panel_top = max(ci_high), .groups = "drop"),
    by = "outcome"
  ) %>%
  mutate(
    x_left = 1,
    x_right = 2,
    x_mid = 1.5,
    y = pmin(0.98, panel_top + 0.06),
    y_tip = y - 0.025,
    label = case_when(
      p_value_holm < 0.001 ~ "***",
      p_value_holm < 0.01 ~ "**",
      TRUE ~ "*"
    )
  )

plot_overall <- plot_overall +
  geom_segment(
    data = overall_brackets,
    aes(x = x_left, xend = x_right, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = overall_brackets,
    aes(x = x_left, xend = x_left, y = y_tip, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = overall_brackets,
    aes(x = x_right, xend = x_right, y = y_tip, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = overall_brackets,
    aes(x = x_mid, y = y + 0.015, label = label),
    inherit.aes = FALSE,
    vjust = 0,
    size = 4
  )

ggsave(
  filename = file.path(output_dir, "overall_accuracy_by_app.png"),
  plot = plot_overall,
  width = 10,
  height = 5,
  dpi = 300
)

plot_skin_tone <- ggplot(skin_tone_summary, aes(x = skin_tone, y = accuracy, color = app, group = app)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.12) +
  facet_wrap(~ outcome, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "DDI accuracy by skin tone group",
    x = "DDI skin tone group",
    y = "Accuracy",
    color = "App"
  )

skin_tone_annotations <- skin_tone_summary %>%
  group_by(outcome, app) %>%
  summarise(y = pmin(0.98, max(ci_high) + 0.08), .groups = "drop") %>%
  left_join(
    skin_tone_ordinal_tests %>%
      transmute(
        app,
        outcome,
        trend_label = if_else(
          !is.na(p_value) & p_value < 0.05,
          paste0("trend ", case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            TRUE ~ "*"
          )),
          NA_character_
        )
      ),
    by = c("app", "outcome")
  ) %>%
  left_join(
    skin_tone_pairwise_tests %>%
      filter(!is.na(p_value_holm), p_value_holm < 0.05) %>%
      mutate(
        pair_piece = paste0(comparison, " ", case_when(
          p_value_holm < 0.001 ~ "***",
          p_value_holm < 0.01 ~ "**",
          TRUE ~ "*"
        ))
      ) %>%
      group_by(app, outcome) %>%
      summarise(
        pair_label = paste(pair_piece, collapse = "\n"),
        .groups = "drop"
      ),
    by = c("app", "outcome")
  ) %>%
  mutate(
    label = case_when(
      !is.na(trend_label) & !is.na(pair_label) ~ paste(trend_label, pair_label, sep = "\n"),
      !is.na(trend_label) ~ trend_label,
      !is.na(pair_label) ~ pair_label,
      TRUE ~ NA_character_
    ),
    x = factor("56", levels = c("12", "34", "56")),
    hjust = if_else(as.character(app) == "ChatGPT", 1.1, -0.1)
  ) %>%
  filter(!is.na(label))

plot_skin_tone <- plot_skin_tone +
  geom_text(
    data = skin_tone_annotations,
    aes(x = x, y = y, label = label, color = app),
    inherit.aes = FALSE,
    hjust = skin_tone_annotations$hjust,
    size = 3.2,
    lineheight = 0.95,
    show.legend = FALSE
  )

ggsave(
  filename = file.path(output_dir, "accuracy_by_skin_tone.png"),
  plot = plot_skin_tone,
  width = 10,
  height = 5,
  dpi = 300
)

plot_disease_by_app <- disease_summary %>%
  mutate(disease_chr = as.character(disease)) %>%
  group_by(disease_chr) %>%
  mutate(max_accuracy_for_order = max(accuracy, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(disease_chr = fct_reorder(disease_chr, max_accuracy_for_order, .desc = TRUE)) %>%
  ggplot(aes(x = accuracy, y = disease_chr, color = app)) +
  geom_point(size = 2.3) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.18) +
  facet_wrap(~ outcome, ncol = 1) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "DDI disease-level accuracy within each app",
    x = "Accuracy",
    y = "Disease",
    color = "App"
  )

ggsave(
  filename = file.path(output_dir, "accuracy_by_disease_within_app.png"),
  plot = plot_disease_by_app,
  width = 11,
  height = 9,
  dpi = 300
)

heatmap_data <- disease_tone_summary %>%
  filter(outcome == "Top-1 diagnosis accuracy") %>%
  mutate(disease = fct_reorder(disease, accuracy, .fun = mean, .desc = TRUE))

plot_heatmap <- ggplot(heatmap_data, aes(x = skin_tone, y = disease, fill = accuracy)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ app, ncol = 2) +
  scale_fill_gradient(low = "#f4efe6", high = "#0f5c5c", labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "DDI top-1 diagnosis accuracy within disease and skin tone",
    x = "DDI skin tone group",
    y = "Disease",
    fill = "Accuracy"
  ) +
  theme(axis.text.y = element_text(size = 8))

ggsave(
  filename = file.path(output_dir, "heatmap_top1_dx_by_disease_skin_tone.png"),
  plot = plot_heatmap,
  width = 11,
  height = 9,
  dpi = 300
)

message("DDI analysis complete.")
message("Outputs written to: ", output_dir)
