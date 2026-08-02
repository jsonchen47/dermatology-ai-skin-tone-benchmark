#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  normalizePath("analysis")
}
repo_dir <- normalizePath(file.path(script_dir, ".."))
data_dir <- file.path(repo_dir, "data", "fitzpatrick")
output_dir <- file.path(repo_dir, "results", "intermediate", "broad_vs_specific")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

regular_path <- file.path(data_dir, "fitzpatrick_analysis_231.csv")
exact_path <- file.path(data_dir, "fitzpatrick_analysis_exact_top1_231.csv")
source_path <- file.path(data_dir, "per_image_application_outputs_231.csv")

app_map <- c(
  a = "AI Skin Scanner",
  m = "Model Dermatol",
  c = "ChatGPT",
  cl = "Claude"
)

app_levels <- unname(app_map)
disease_levels <- c("basal cell carcinoma", "benign nevus", "melanoma spectrum", "squamous cell carcinoma")

app_palette <- c(
  "AI Skin Scanner" = "#D55E00",
  "Model Dermatol" = "#0072B2",
  "ChatGPT" = "#009E73",
  "Claude" = "#CC79A7"
)

mapping_palette <- c(
  "Broad allowed" = "#9CC2E5",
  "Specific only" = "#1F4E79"
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

disease_group <- function(label) {
  case_when(
    label %in% c("melanoma", "malignant melanoma") ~ "melanoma spectrum",
    label %in% c("congenital nevus", "halo nevus") ~ "benign nevus",
    label == "basal cell carcinoma" ~ "basal cell carcinoma",
    label == "squamous cell carcinoma" ~ "squamous cell carcinoma",
    TRUE ~ NA_character_
  )
}

regular <- read_csv(regular_path, show_col_types = FALSE) %>%
  mutate(disease_group = disease_group(label)) %>%
  filter(!is.na(disease_group))

exact <- read_csv(exact_path, show_col_types = FALSE) %>%
  mutate(disease_group = disease_group(label)) %>%
  filter(!is.na(disease_group))

source <- read_csv(source_path, show_col_types = FALSE) %>%
  mutate(disease_group = disease_group(label)) %>%
  filter(!is.na(disease_group)) %>%
  select(md5hash, disease_group, a_dx_1, m_dx_1, c_dx_1, cl_dx_1)

make_long <- function(df, suffix, mapping_name) {
  df %>%
    transmute(
      md5hash,
      disease_group,
      `AI Skin Scanner` = as.integer(.data[[paste0("a_", suffix)]]),
      `Model Dermatol` = as.integer(.data[[paste0("m_", suffix)]]),
      `ChatGPT` = as.integer(.data[[paste0("c_", suffix)]]),
      `Claude` = as.integer(.data[[paste0("cl_", suffix)]])
    ) %>%
    pivot_longer(cols = all_of(app_levels), names_to = "app", values_to = "correct") %>%
    mutate(mapping = mapping_name)
}

regular_long <- make_long(regular, "dx_1_correct", "Broad allowed")
exact_long <- make_long(exact, "dx_1_correct", "Specific only")

combined <- bind_rows(regular_long, exact_long) %>%
  mutate(
    disease_group = factor(disease_group, levels = disease_levels),
    app = factor(app, levels = app_levels),
    mapping = factor(mapping, levels = c("Broad allowed", "Specific only"))
  )

summary_tbl <- combined %>%
  group_by(disease_group, app, mapping) %>%
  summarise(
    n = n(),
    correct_n = sum(correct),
    accuracy = correct_n / n,
    ci_low = wilson_ci(correct_n, n)[1],
    ci_high = wilson_ci(correct_n, n)[2],
    .groups = "drop"
  )

delta_tbl <- summary_tbl %>%
  select(disease_group, app, mapping, accuracy, correct_n, n) %>%
  pivot_wider(names_from = mapping, values_from = c(accuracy, correct_n, n)) %>%
  mutate(
    drop_pts = 100 * (`accuracy_Specific only` - `accuracy_Broad allowed`),
    drop_n = `correct_n_Broad allowed` - `correct_n_Specific only`,
    label = sprintf("%+.1f pts", drop_pts)
  )

write_csv(summary_tbl, file.path(output_dir, "broad_vs_specific_top1_by_disease_app.csv"))
write_csv(delta_tbl, file.path(output_dir, "broad_vs_specific_top1_deltas_by_disease_app.csv"))

regular_exact_join <- regular %>%
  select(md5hash, disease_group, starts_with("a_dx_1_correct"), starts_with("m_dx_1_correct"), starts_with("c_dx_1_correct"), starts_with("cl_dx_1_correct")) %>%
  rename_with(~ str_replace(.x, "_dx_1_correct$", "_regular"), ends_with("dx_1_correct")) %>%
  left_join(
    exact %>%
      select(md5hash, starts_with("a_dx_1_correct"), starts_with("m_dx_1_correct"), starts_with("c_dx_1_correct"), starts_with("cl_dx_1_correct")) %>%
      rename_with(~ str_replace(.x, "_dx_1_correct$", "_exact"), ends_with("dx_1_correct")),
    by = "md5hash"
  ) %>%
  left_join(source, by = c("md5hash", "disease_group"))

removed_terms <- bind_rows(lapply(names(app_map), function(prefix) {
  regular_col <- paste0(prefix, "_regular")
  exact_col <- paste0(prefix, "_exact")
  dx_col <- paste0(prefix, "_dx_1")
  regular_exact_join %>%
    filter(.data[[regular_col]] == 1, .data[[exact_col]] == 0) %>%
    transmute(
      disease_group,
      app = app_map[[prefix]],
      removed_top1_term = .data[[dx_col]]
    )
})) %>%
  mutate(
    disease_group = factor(disease_group, levels = disease_levels),
    app = factor(app, levels = app_levels)
  )

removed_terms_summary <- removed_terms %>%
  count(disease_group, app, removed_top1_term, sort = TRUE, name = "n_removed") %>%
  group_by(disease_group, app) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 5) %>%
  ungroup()

removed_terms_labels <- removed_terms_summary %>%
  group_by(disease_group, app) %>%
  summarise(
    removed_terms = paste0(removed_top1_term, " (", n_removed, ")", collapse = "; "),
    .groups = "drop"
  )

write_csv(removed_terms_summary, file.path(output_dir, "removed_broad_terms_by_disease_app.csv"))
write_csv(removed_terms_labels, file.path(output_dir, "removed_broad_terms_labels_by_disease_app.csv"))

p1 <- ggplot(summary_tbl, aes(x = app, y = accuracy, fill = mapping)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    position = position_dodge(width = 0.72),
    width = 0.14
  ) +
  geom_text(
    data = delta_tbl,
    aes(x = app, y = 0.93, label = label),
    inherit.aes = FALSE,
    size = 3.1,
    fontface = "bold",
    color = "#374151"
  ) +
  facet_wrap(~ disease_group, nrow = 1) +
  scale_fill_manual(values = mapping_palette) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.98)) +
  labs(
    title = "Counting broad labels as correct can materially change apparent top-1 performance",
    subtitle = "Per disease and per app: Broad allowed = regular mapping; Specific only = broad terms like skin cancer removed. Labels show specific-only minus broad-allowed change.",
    x = NULL,
    y = "Top-1 diagnosis accuracy",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "fig_broad_vs_specific_top1_by_disease_app.png"), p1, width = 15.2, height = 6.0, dpi = 300)

heatmap_tbl <- delta_tbl %>%
  mutate(
    cell_label = paste0(label, "\n", drop_n, " cases"),
    disease_group = factor(disease_group, levels = disease_levels),
    app = factor(app, levels = app_levels)
  )

p2 <- ggplot(heatmap_tbl, aes(x = app, y = disease_group, fill = drop_pts)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = cell_label), size = 3.0, lineheight = 0.95) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "#F7F7F7",
    high = "#2166AC",
    midpoint = 0,
    labels = function(x) sprintf("%+.0f pts", x)
  ) +
  labs(
    title = "The effect of removing broad terms is concentrated in specific disease-app combinations",
    subtitle = "Negative values mean top-1 accuracy fell when broad labels were no longer counted as correct",
    x = NULL,
    y = NULL,
    fill = "Change"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "fig_broad_vs_specific_drop_heatmap.png"), p2, width = 10.8, height = 5.8, dpi = 300)
