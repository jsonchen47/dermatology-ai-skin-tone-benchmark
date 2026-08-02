#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  normalizePath("analysis")
}
repo_dir <- normalizePath(file.path(script_dir, ".."))
csv_path <- file.path(repo_dir, "data", "fitzpatrick", "cohort_manifest_231.csv")
output_dir <- file.path(repo_dir, "results", "power_context")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(csv_path)) {
  stop("Could not find CSV at: ", csv_path)
}

fitz <- read_csv(csv_path, show_col_types = FALSE) %>%
  rename(label = source_label)

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

find_mde_equal <- function(n, p1 = 0.5, power = 0.8, sig.level = 0.05) {
  vals <- c(
    seq(max(0.001, p1 - 0.6), p1 - 0.001, by = 0.001),
    seq(p1 + 0.001, min(0.999, p1 + 0.6), by = 0.001)
  )
  pw <- vapply(
    vals,
    function(p2) power.prop.test(n = n, p1 = p1, p2 = p2, sig.level = sig.level)$power,
    numeric(1)
  )
  ok <- which(pw >= power)
  if (length(ok) == 0) {
    return(c(detectable_p2 = NA_real_, detectable_abs_difference = NA_real_))
  }
  best <- ok[which.min(abs(vals[ok] - p1))]
  c(
    detectable_p2 = vals[best],
    detectable_abs_difference = abs(vals[best] - p1)
  )
}

subset_counts <- fitz %>%
  mutate(
    fitzpatrick_scale = suppressWarnings(as.integer(fitzpatrick_scale)),
    disease_group = disease_group(label),
    grouped_tone = tone_group(fitzpatrick_scale)
  ) %>%
  filter(!is.na(disease_group), !is.na(grouped_tone), !is.na(fitzpatrick_scale))

grouped_counts <- subset_counts %>%
  count(disease_group, grouped_tone, name = "n") %>%
  tidyr::pivot_wider(names_from = grouped_tone, values_from = n, values_fill = 0) %>%
  mutate(total = `I-II` + `III-IV` + `V-VI`) %>%
  arrange(factor(disease_group, levels = c(
    "basal cell carcinoma",
    "squamous cell carcinoma",
    "melanoma spectrum",
    "benign nevus"
  )))

level_counts <- subset_counts %>%
  count(disease_group, fitzpatrick_scale, name = "n") %>%
  tidyr::pivot_wider(names_from = fitzpatrick_scale, values_from = n, values_fill = 0) %>%
  mutate(total = `1` + `2` + `3` + `4` + `5` + `6`) %>%
  arrange(factor(disease_group, levels = c(
    "basal cell carcinoma",
    "squamous cell carcinoma",
    "melanoma spectrum",
    "benign nevus"
  )))

limiting_grouped <- grouped_counts %>%
  pivot_longer(cols = c(`I-II`, `III-IV`, `V-VI`), names_to = "grouped_tone", values_to = "n") %>%
  arrange(n, disease_group, grouped_tone)

limiting_levels <- level_counts %>%
  pivot_longer(cols = c(`1`, `2`, `3`, `4`, `5`, `6`), names_to = "fitzpatrick_scale", values_to = "n") %>%
  arrange(n, disease_group, fitzpatrick_scale)

power_table <- bind_rows(
  tibble(
    scenario = "Overall grouped tones (~70 per group)",
    n_per_group = 70,
    baseline_accuracy = c(0.30, 0.50, 0.70)
  ),
  tibble(
    scenario = "Disease-specific grouped tones (~20 per group)",
    n_per_group = 20,
    baseline_accuracy = c(0.30, 0.50, 0.70)
  ),
  tibble(
    scenario = "Small 6-level cells (~10 per group)",
    n_per_group = 10,
    baseline_accuracy = c(0.30, 0.50, 0.70)
  )
) %>%
  mutate(calc = map2(n_per_group, baseline_accuracy, ~ find_mde_equal(.x, .y))) %>%
  mutate(
    detectable_p2 = map_dbl(calc, ~ .x[["detectable_p2"]]),
    detectable_abs_difference = map_dbl(calc, ~ .x[["detectable_abs_difference"]])
  ) %>%
  mutate(
    baseline_accuracy_pct = round(baseline_accuracy * 100, 1),
    detectable_p2_pct = round(detectable_p2 * 100, 1),
    detectable_abs_difference_pct_points = round(detectable_abs_difference * 100, 1)
  ) %>%
  select(
    scenario,
    n_per_group,
    baseline_accuracy_pct,
    detectable_p2_pct,
    detectable_abs_difference_pct_points
  )

write_csv(grouped_counts, file.path(output_dir, "fitzpatrick_grouped_tone_counts_for_power_context.csv"))
write_csv(level_counts, file.path(output_dir, "fitzpatrick_6level_tone_counts_for_power_context.csv"))
write_csv(limiting_grouped, file.path(output_dir, "fitzpatrick_limiting_grouped_tone_cells.csv"))
write_csv(limiting_levels, file.path(output_dir, "fitzpatrick_limiting_6level_tone_cells.csv"))
write_csv(power_table, file.path(output_dir, "quick_power_analysis_table.csv"))

cat("Wrote grouped counts, level counts, limiting cells, and quick power table to:\n")
cat(output_dir, "\n")
