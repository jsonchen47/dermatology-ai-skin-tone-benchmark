#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(purrr)
  library(stringr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  normalizePath("analysis")
}
repo_dir <- normalizePath(file.path(script_dir, ".."))
output_dir <- file.path(repo_dir, "results", "intermediate", "full_231_significance")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fitz_path <- file.path(repo_dir, "data", "fitzpatrick", "fitzpatrick_analysis_231.csv")

app_map <- c(
  a = "AI Skin Scanner",
  m = "Model Dermatol",
  c = "ChatGPT",
  cl = "Claude"
)

app_levels <- unname(app_map)
metric_levels <- c("Top-1 diagnosis", "Top-3 diagnosis", "Benign/malignant triage")
disease_levels <- c("basal cell carcinoma", "benign nevus", "melanoma spectrum", "squamous cell carcinoma")
tone_levels <- c("I-II", "III-IV", "V-VI")

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

wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf.level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2))) / denom
  c(max(0, center - half), min(1, center + half))
}

mean_ci <- function(x) {
  n <- sum(!is.na(x))
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(n)
  lo <- m - qt(0.975, df = n - 1) * se
  hi <- m + qt(0.975, df = n - 1) * se
  c(mean = m, ci_low = max(0, lo), ci_high = min(1, hi), n = n)
}

fmt_p <- function(p) {
  ifelse(
    p < 0.001,
    formatC(p, format = "e", digits = 1),
    sub("^0", "", sprintf("%.3f", p))
  )
}

disease_group <- function(label) {
  case_when(
    label %in% c("melanoma", "malignant melanoma") ~ "melanoma spectrum",
    label %in% c("congenital nevus", "halo nevus") ~ "benign nevus",
    label == "basal cell carcinoma" ~ "basal cell carcinoma",
    label == "squamous cell carcinoma" ~ "squamous cell carcinoma",
    TRUE ~ NA_character_
  )
}

tone_group <- function(scale) {
  case_when(
    scale %in% c(1, 2) ~ "I-II",
    scale %in% c(3, 4) ~ "III-IV",
    scale %in% c(5, 6) ~ "V-VI",
    TRUE ~ NA_character_
  )
}

cochran_q_p <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  k <- ncol(mat)
  col_success <- colSums(mat)
  row_success <- rowSums(mat)
  total_success <- sum(col_success)
  q_stat <- ((k - 1) * (k * sum(col_success^2) - total_success^2)) /
    ((k * total_success) - sum(row_success^2))
  pchisq(q_stat, df = k - 1, lower.tail = FALSE)
}

pairwise_mcnemar <- function(df, suffix) {
  combn(names(app_map), 2, simplify = FALSE, FUN = function(pair) {
    x <- as.integer(df[[paste0(pair[1], "_", suffix)]])
    y <- as.integer(df[[paste0(pair[2], "_", suffix)]])
    keep <- !is.na(x) & !is.na(y)
    x <- x[keep]
    y <- y[keep]
    d10 <- sum(x == 1 & y == 0)
    d01 <- sum(x == 0 & y == 1)
    p <- mcnemar.test(matrix(c(0, d10, d01, 0), nrow = 2), correct = TRUE)$p.value
    tibble(
      app_1 = app_map[[pair[1]]],
      app_2 = app_map[[pair[2]]],
      discordant_10 = d10,
      discordant_01 = d01,
      p_value = p,
      better_app = case_when(
        d10 > d01 ~ app_map[[pair[1]]],
        d01 > d10 ~ app_map[[pair[2]]],
        TRUE ~ "Tie"
      )
    )
  }) %>% bind_rows()
}

fitz <- read_csv(fitz_path, show_col_types = FALSE) %>%
  mutate(
    fitzpatrick_scale = as.integer(fitzpatrick_scale),
    disease_group = disease_group(label),
    tone_group = tone_group(fitzpatrick_scale)
  ) %>%
  filter(!is.na(disease_group), !is.na(tone_group))

make_long <- function(df, suffix, metric_name) {
  df %>%
    transmute(
      disease_group,
      tone_group,
      `AI Skin Scanner` = as.integer(.data[[paste0("a_", suffix)]]),
      `Model Dermatol` = as.integer(.data[[paste0("m_", suffix)]]),
      `ChatGPT` = as.integer(.data[[paste0("c_", suffix)]]),
      `Claude` = as.integer(.data[[paste0("cl_", suffix)]])
    ) %>%
    pivot_longer(cols = all_of(app_levels), names_to = "app", values_to = "correct") %>%
    mutate(
      metric = metric_name,
      app_family = if_else(app %in% c("ChatGPT", "Claude"), "LLMs", "Smartphone apps")
    )
}

top1_long <- make_long(fitz, "dx_1_correct", "Top-1 diagnosis")
top3_long <- make_long(fitz, "any_dx_correct", "Top-3 diagnosis")
mal_long <- make_long(fitz, "dx_1_malignancy_correct", "Benign/malignant triage")
all_long <- bind_rows(top1_long, top3_long, mal_long) %>%
  mutate(
    metric = factor(metric, levels = metric_levels),
    app = factor(app, levels = app_levels),
    disease_group = factor(disease_group, levels = disease_levels),
    tone_group = factor(tone_group, levels = tone_levels)
  )

overall_by_app <- all_long %>%
  group_by(metric, app) %>%
  summarise(
    n = n(),
    correct_n = sum(correct),
    accuracy = correct_n / n,
    ci_low = wilson_ci(correct_n, n)[1],
    ci_high = wilson_ci(correct_n, n)[2],
    .groups = "drop"
  ) %>%
  mutate(metric = factor(metric, levels = metric_levels))

overall_tests <- bind_rows(
  pairwise_mcnemar(fitz, "dx_1_correct") %>% mutate(metric = "Top-1 diagnosis"),
  pairwise_mcnemar(fitz, "any_dx_correct") %>% mutate(metric = "Top-3 diagnosis"),
  pairwise_mcnemar(fitz, "dx_1_malignancy_correct") %>% mutate(metric = "Benign/malignant triage")
) %>%
  group_by(metric) %>%
  mutate(p_value_holm = p.adjust(p_value, method = "holm")) %>%
  ungroup()

overall_omnibus <- tibble(
  metric = metric_levels,
  cochran_q_p = c(
    cochran_q_p(fitz[, c("a_dx_1_correct", "m_dx_1_correct", "c_dx_1_correct", "cl_dx_1_correct")]),
    cochran_q_p(fitz[, c("a_any_dx_correct", "m_any_dx_correct", "c_any_dx_correct", "cl_any_dx_correct")]),
    cochran_q_p(fitz[, c("a_dx_1_malignancy_correct", "m_dx_1_malignancy_correct", "c_dx_1_malignancy_correct", "cl_dx_1_malignancy_correct")])
  )
)

overall_ann <- overall_tests %>%
  filter(p_value_holm < 0.05) %>%
  mutate(
    worse_app = if_else(better_app == app_1, app_2, app_1),
    line = paste0(better_app, " > ", worse_app, ", Holm P=", fmt_p(p_value_holm))
  ) %>%
  group_by(metric) %>%
  summarise(sig_lines = paste(line, collapse = "\n"), .groups = "drop") %>%
  right_join(overall_omnibus, by = "metric") %>%
  mutate(
    metric = factor(metric, levels = metric_levels),
    sig_lines = replace_na(sig_lines, "No significant pairwise differences"),
    label = case_when(
      metric == "Top-1 diagnosis" ~ paste0(
        "Overall P<.001",
        "\nModel Dermatol > all others",
        "\nHolm P=.003 or lower"
      ),
      metric == "Top-3 diagnosis" ~ paste0(
        "Overall P<.001",
        "\nChatGPT > AI Skin Scanner, P<.001",
        "\nChatGPT > Claude, P=.038"
      ),
      metric == "Benign/malignant triage" ~ paste0(
        "Overall P<.001",
        "\nClaude > all others",
        "\nHolm P=.039 or lower"
      )
    ),
    y = 1.08
  )

write_csv(overall_by_app, file.path(output_dir, "overall_by_app_with_significance_full231.csv"))
write_csv(overall_tests, file.path(output_dir, "overall_pairwise_tests_with_significance_full231.csv"))
write_csv(overall_omnibus, file.path(output_dir, "overall_omnibus_tests_full231.csv"))

family_metric_info <- list(
  "Top-1 diagnosis" = c("a_dx_1_correct", "m_dx_1_correct", "c_dx_1_correct", "cl_dx_1_correct"),
  "Top-3 diagnosis" = c("a_any_dx_correct", "m_any_dx_correct", "c_any_dx_correct", "cl_any_dx_correct"),
  "Benign/malignant triage" = c("a_dx_1_malignancy_correct", "m_dx_1_malignancy_correct", "c_dx_1_malignancy_correct", "cl_dx_1_malignancy_correct")
)

family_tests_raw <- bind_rows(lapply(names(family_metric_info), function(metric_name) {
  cols <- family_metric_info[[metric_name]]
  phone <- (fitz[[cols[1]]] + fitz[[cols[2]]]) / 2
  llm <- (fitz[[cols[3]]] + fitz[[cols[4]]]) / 2
  wt <- wilcox.test(llm, phone, paired = TRUE, exact = FALSE)
  llm_ci <- mean_ci(llm)
  phone_ci <- mean_ci(phone)
  bind_rows(
    tibble(metric = metric_name, family = "LLMs", mean = llm_ci["mean"], ci_low = llm_ci["ci_low"], ci_high = llm_ci["ci_high"]),
    tibble(metric = metric_name, family = "Smartphone apps", mean = phone_ci["mean"], ci_low = phone_ci["ci_low"], ci_high = phone_ci["ci_high"])
  ) %>%
    mutate(
      test_p = wt$p.value,
      llm_mean = mean(llm),
      phone_mean = mean(phone)
    )
})) %>%
  mutate(metric = factor(metric, levels = metric_levels))

family_tests <- family_tests_raw %>%
  distinct(metric, test_p, llm_mean, phone_mean) %>%
  mutate(
    p_holm = p.adjust(test_p, method = "holm"),
    winner = if_else(llm_mean > phone_mean, "LLMs", "Smartphone apps"),
    sig_label = case_when(
      p_holm < 0.05 ~ paste0(winner, " higher, Holm P=", fmt_p(p_holm)),
      TRUE ~ paste0("Not significant, Holm P=", fmt_p(p_holm))
    )
  )

family_plot <- family_tests_raw %>%
  select(metric, family, mean, ci_low, ci_high) %>%
  mutate(
    family = factor(family, levels = c("Smartphone apps", "LLMs")),
    metric = factor(metric, levels = metric_levels)
  )

write_csv(family_plot, file.path(output_dir, "family_metric_plot_values_full231.csv"))
write_csv(family_tests, file.path(output_dir, "family_metric_significance_full231.csv"))

disease_top1_by_app <- top1_long %>%
  group_by(disease_group, app) %>%
  summarise(
    n = n(),
    correct_n = sum(correct),
    accuracy = correct_n / n,
    ci_low = wilson_ci(correct_n, n)[1],
    ci_high = wilson_ci(correct_n, n)[2],
    .groups = "drop"
  ) %>%
  mutate(
    disease_group = factor(disease_group, levels = disease_levels),
    app = factor(app, levels = app_levels)
  )

disease_pairwise <- map_dfr(disease_levels, function(disease_name) {
  df_d <- fitz %>% filter(disease_group == disease_name)
  pairwise_mcnemar(df_d, "dx_1_correct") %>%
    mutate(disease_group = disease_name)
}) %>%
  group_by(disease_group) %>%
  mutate(p_value_holm = p.adjust(p_value, method = "holm")) %>%
  ungroup() %>%
  mutate(
    disease_group = factor(disease_group, levels = disease_levels)
  )

disease_sig_labels <- disease_pairwise %>%
  filter(p_value_holm < 0.05) %>%
  mutate(line = paste0(better_app, " > ", if_else(better_app == app_1, app_2, app_1), ", P=", fmt_p(p_value_holm))) %>%
  group_by(disease_group) %>%
  summarise(sig_lines = paste(line, collapse = "\n"), .groups = "drop")

disease_winner_labels <- disease_top1_by_app %>%
  group_by(disease_group) %>%
  slice_max(order_by = accuracy, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(disease_group, top_app = as.character(app))

disease_summary_labels <- disease_winner_labels %>%
  left_join(disease_sig_labels, by = "disease_group") %>%
  mutate(
    label = case_when(
      disease_group == "basal cell carcinoma" ~ paste0(
        "Model Dermatol > ChatGPT and Claude\n",
        "AI Skin Scanner and ChatGPT > Claude\n",
        "all Holm P=.018 or lower"
      ),
      disease_group == "benign nevus" ~ paste0(
        "AI Skin Scanner lower than all others\n",
        "all Holm P=.023 or lower"
      ),
      disease_group == "melanoma spectrum" ~ paste0(
        "No significant pairwise differences\n",
        "after Holm correction"
      ),
      disease_group == "squamous cell carcinoma" ~ paste0(
        "Model Dermatol and Claude\n",
        "> AI Skin Scanner\n",
        "both Holm P<.001"
      ),
      TRUE ~ paste0("Highest raw accuracy: ", top_app)
    )
  )

write_csv(disease_top1_by_app, file.path(output_dir, "disease_top1_by_app_full231.csv"))
write_csv(disease_pairwise, file.path(output_dir, "disease_pairwise_top1_full231.csv"))
write_csv(disease_summary_labels, file.path(output_dir, "disease_summary_labels_full231.csv"))

tone_top1 <- top1_long %>%
  group_by(app, tone_group) %>%
  summarise(
    n = n(),
    correct_n = sum(correct),
    accuracy = correct_n / n,
    ci_low = wilson_ci(correct_n, n)[1],
    ci_high = wilson_ci(correct_n, n)[2],
    .groups = "drop"
  ) %>%
  mutate(
    app = factor(app, levels = app_levels),
    tone_group = factor(tone_group, levels = tone_levels)
  )

pooled_tone_tests <- map_dfr(names(app_map), function(prefix) {
  app_name <- app_map[[prefix]]
  tmp <- fitz %>%
    count(tone_group, correct = .data[[paste0(prefix, "_dx_1_correct")]]) %>%
    pivot_wider(names_from = correct, values_from = n, values_fill = 0) %>%
    arrange(match(tone_group, tone_levels))

  tone_mat <- as.matrix(tmp[, c("1", "0")])
  rownames(tone_mat) <- tmp$tone_group

  tibble(
    app = app_name,
    grouped_overall_p = suppressWarnings(chisq.test(tone_mat)$p.value),
    grouped_ordinal_p = prop.trend.test(as.numeric(tone_mat[, "1"]), rowSums(tone_mat), score = c(1, 2, 3))$p.value
  )
}) %>%
  mutate(
    grouped_overall_p_holm = p.adjust(grouped_overall_p, method = "holm"),
    grouped_ordinal_p_holm = p.adjust(grouped_ordinal_p, method = "holm")
  )

claude_tone_label <- pooled_tone_tests %>%
  filter(app == "Claude") %>%
  mutate(
    label = paste0(
      "Significant\n",
      "Overall P=", fmt_p(grouped_overall_p_holm), "\n",
      "Ordinal P=", fmt_p(grouped_ordinal_p_holm)
    )
  )

write_csv(tone_top1, file.path(output_dir, "tone_top1_by_app_full231.csv"))
write_csv(pooled_tone_tests, file.path(output_dir, "pooled_tone_tests_top1_full231.csv"))

p1 <- ggplot(overall_by_app, aes(x = app, y = accuracy, fill = app)) +
  geom_col(width = 0.68) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
  facet_wrap(~ metric, nrow = 1) +
  geom_text(
    data = overall_ann,
    aes(x = 2.5, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.8,
    lineheight = 1.0,
    color = "#374151"
  ) +
  scale_fill_manual(values = app_palette) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.17),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "Individual app performance differs significantly,\nand the ranking changes by endpoint",
    subtitle = "Full 231-image cohort; annotations summarize omnibus and Holm-corrected pairwise results",
    x = NULL,
    y = "Accuracy"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", lineheight = 1.0),
    plot.subtitle = element_text(size = 11.5),
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 12),
    plot.margin = margin(10, 10, 10, 10)
  )
ggsave(file.path(output_dir, "fig1_individual_apps_by_metric_with_significance.png"), p1, width = 9.2, height = 6.5, dpi = 600)
ggsave(file.path(output_dir, "fig1_individual_apps_by_metric_with_significance.pdf"), p1, width = 9.2, height = 6.5)

family_plot_display <- family_plot %>%
  mutate(
    metric_plot = factor(as.character(metric), levels = rev(metric_levels)),
    value_label = percent(mean, accuracy = 0.1)
  )

family_segments <- family_plot_display %>%
  select(metric_plot, family, mean) %>%
  pivot_wider(names_from = family, values_from = mean)

family_test_display <- family_tests %>%
  mutate(metric_plot = factor(as.character(metric), levels = rev(metric_levels)))

p2 <- ggplot() +
  geom_segment(
    data = family_segments,
    aes(
      x = `Smartphone apps`,
      xend = LLMs,
      y = metric_plot,
      yend = metric_plot
    ),
    linewidth = 1.8,
    color = "#D1D5DB"
  ) +
  geom_errorbar(
    data = family_plot_display,
    aes(
      x = mean,
      y = metric_plot,
      xmin = ci_low,
      xmax = ci_high,
      color = family
    ),
    orientation = "y",
    width = 0.16,
    linewidth = 0.75
  ) +
  geom_point(
    data = family_plot_display,
    aes(x = mean, y = metric_plot, color = family),
    size = 5
  ) +
  geom_text(
    data = family_plot_display,
    aes(x = mean, y = metric_plot, label = value_label),
    nudge_y = -0.24,
    size = 3.8,
    color = "#111827"
  ) +
  geom_text(
    data = family_test_display,
    aes(x = 1.01, y = metric_plot, label = sig_label),
    hjust = 1,
    vjust = -1.15,
    size = 3.8,
    color = "#374151"
  ) +
  scale_color_manual(values = family_palette) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0.3, 1.03),
    breaks = seq(0.4, 1, 0.2)
  ) +
  labs(
    title = "LLM-vs-smartphone differences are significant only\nfor benign/malignant triage",
    subtitle = paste0(
      "Family score = mean correctness of the two apps in each family for each image;\n",
      "paired Wilcoxon P values are Holm-corrected across the 3 endpoints"
    ),
    x = "Family accuracy (95% CI)",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = 16, face = "bold", lineheight = 1.0),
    plot.subtitle = element_text(size = 10.5),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )
ggsave(file.path(output_dir, "fig2_family_comparison_by_metric_with_significance.png"), p2, width = 9.2, height = 5.8, dpi = 600)
ggsave(file.path(output_dir, "fig2_family_comparison_by_metric_with_significance_doc.png"), p2, width = 9.2, height = 6.1, dpi = 300)
ggsave(file.path(output_dir, "fig2_family_comparison_by_metric_with_significance.pdf"), p2, width = 9.2, height = 5.8)

p3 <- ggplot(disease_top1_by_app, aes(x = app, y = accuracy, fill = app)) +
  geom_col(width = 0.68) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15) +
  facet_wrap(~ disease_group, ncol = 2) +
  geom_text(
    data = disease_summary_labels,
    aes(x = 2.5, y = 1.07, label = label),
    inherit.aes = FALSE,
    size = 4.0,
    lineheight = 1.0,
    color = "#374151"
  ) +
  scale_fill_manual(values = app_palette) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.16),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "Disease-specific performance patterns are large,\nbut no application dominates every disease",
    subtitle = "Regular top-1 diagnosis accuracy; annotations summarize Holm-corrected within-disease comparisons",
    x = NULL,
    y = "Top-1 diagnosis accuracy"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 16, face = "bold", lineheight = 1.0),
    plot.subtitle = element_text(size = 10.5),
    axis.text.x = element_text(angle = 20, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    strip.text = element_text(face = "bold", size = 11),
    plot.margin = margin(10, 10, 10, 10)
  )
ggsave(file.path(output_dir, "fig3_disease_performance_with_significance.png"), p3, width = 9.2, height = 9.0, dpi = 600)
ggsave(file.path(output_dir, "fig3_disease_performance_with_significance_doc.png"), p3, width = 9.2, height = 9.0, dpi = 300)
ggsave(file.path(output_dir, "fig3_disease_performance_with_significance.pdf"), p3, width = 9.2, height = 9.0)

tone_other <- tone_top1 %>% filter(app != "Claude")
tone_claude <- tone_top1 %>% filter(app == "Claude")

p4 <- ggplot(tone_top1, aes(x = tone_group, y = accuracy, group = app, color = app)) +
  geom_line(
    data = tone_other,
    linewidth = 0.9,
    alpha = 0.42
  ) +
  geom_errorbar(
    data = tone_other,
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.07,
    linewidth = 0.55,
    alpha = 0.42
  ) +
  geom_point(
    data = tone_other,
    size = 3.4,
    alpha = 0.62
  ) +
  geom_line(
    data = tone_claude,
    linewidth = 1.7
  ) +
  geom_errorbar(
    data = tone_claude,
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.08,
    linewidth = 0.85
  ) +
  geom_point(
    data = tone_claude,
    size = 5
  ) +
  geom_text(
    data = tone_claude,
    aes(
      y = pmin(ci_high + 0.035, 0.89),
      label = percent(accuracy, accuracy = 0.1)
    ),
    size = 3.8,
    color = "#111827"
  ) +
  geom_text(
    data = claude_tone_label,
    aes(
      x = "V-VI",
      y = 0.98,
      label = paste0(
        "Claude: overall Holm P=",
        fmt_p(grouped_overall_p_holm),
        "; ordinal Holm P=",
        fmt_p(grouped_ordinal_p_holm)
      )
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = 4.0,
    color = "#374151"
  ) +
  scale_color_manual(values = app_palette) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.04),
    breaks = seq(0, 1, 0.25)
  ) +
  labs(
    title = "Only Claude shows a significant pooled decline\nacross grouped skin-tone strata",
    subtitle = "Regular top-1 diagnosis accuracy by grouped Fitzpatrick tone with 95% confidence intervals",
    x = "Grouped Fitzpatrick skin tone",
    y = "Top-1 diagnosis accuracy",
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(size = 16, face = "bold", lineheight = 1.0),
    plot.subtitle = element_text(size = 10.5),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )
ggsave(file.path(output_dir, "fig4_skin_tone_group_performance_with_claude_significance.png"), p4, width = 9.2, height = 5.8, dpi = 600)
ggsave(file.path(output_dir, "fig4_skin_tone_group_performance_with_claude_significance_doc.png"), p4, width = 9.2, height = 5.0, dpi = 300)
ggsave(file.path(output_dir, "fig4_skin_tone_group_performance_with_claude_significance.pdf"), p4, width = 9.2, height = 5.8)
