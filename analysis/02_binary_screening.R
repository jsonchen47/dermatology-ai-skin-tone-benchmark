#!/usr/bin/env Rscript

# Binary screening analysis for Fitzpatrick dataset:
# - malignant vs benign
# - melanoma vs rest
# - basal cell carcinoma vs rest
# - squamous cell carcinoma vs rest
# Produces pooled and tone-stratified confusion metrics, AUROC/AUPRC within app,
# and inferential summaries for between-app and within-app subgroup differences.

required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "purrr", "ggplot2", "scales"
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
  library(purrr)
  library(ggplot2)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  normalizePath("analysis")
}
repo_dir <- normalizePath(file.path(script_dir, ".."))
data_dir <- file.path(repo_dir, "data", "fitzpatrick")

results_path <- file.path(data_dir, "per_image_application_outputs_231.csv")
truth_diag_path <- file.path(data_dir, "truth_correct_diagnosis.csv")
truth_benmal_path <- file.path(data_dir, "truth_malignant_or_benign.csv")
output_dir <- file.path(repo_dir, "results", "binary_screening")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

apps <- c(
  a = "AI Skin Scanner",
  m = "Model Dermatol",
  c = "ChatGPT",
  cl = "Claude"
)

app_levels <- unname(apps)

targets <- c("malignant_vs_benign", "melanoma_vs_rest", "basal_cell_carcinoma_vs_rest", "squamous_cell_carcinoma_vs_rest")

target_labels <- c(
  malignant_vs_benign = "Malignant vs Benign",
  melanoma_vs_rest = "Melanoma vs Rest",
  basal_cell_carcinoma_vs_rest = "Basal Cell Carcinoma vs Rest",
  squamous_cell_carcinoma_vs_rest = "Squamous Cell Carcinoma vs Rest"
)

metric_labels <- c(
  sensitivity = "Sensitivity",
  specificity = "Specificity",
  ppv = "PPV",
  npv = "NPV",
  auroc = "AUROC",
  auprc = "AUPRC"
)

app_palette <- c(
  "AI Skin Scanner" = "#D55E00",
  "Model Dermatol" = "#0072B2",
  "ChatGPT" = "#009E73",
  "Claude" = "#CC79A7"
)

fitz_group_3 <- function(x) {
  case_when(
    x %in% c("1", "2") ~ "I-II",
    x %in% c("3", "4") ~ "III-IV",
    x %in% c("5", "6") ~ "V-VI",
    TRUE ~ "Unknown"
  )
}

group_truth_disease <- function(x) {
  case_when(
    x %in% c("congenital nevus", "halo nevus") ~ "benign nevus",
    x %in% c("malignant melanoma", "melanoma", "lentigo maligna") ~ "melanoma spectrum",
    TRUE ~ x
  )
}

normalize_benmal <- function(x) {
  case_when(
    str_to_lower(str_squish(x)) %in% c("b", "benign") ~ "benign",
    str_to_lower(str_squish(x)) %in% c("m", "malignant") ~ "malignant",
    TRUE ~ str_to_lower(str_squish(x))
  )
}

split_terms <- function(x) {
  x %>%
    str_split("\\|") %>%
    purrr::pluck(1) %>%
    str_squish() %>%
    purrr::discard(~ .x == "")
}

parse_score <- function(app_code, value) {
  value <- str_squish(as.character(value))
  if (is.na(value) || value == "" || toupper(value) == "N/A") {
    return(NA_real_)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed)) {
    return(NA_real_)
  }

  norm <- case_when(
    app_code == "a" ~ parsed,
    app_code == "m" & parsed > 1 ~ parsed / 100,
    app_code == "m" ~ parsed,
    app_code %in% c("c", "cl") ~ parsed / 100,
    TRUE ~ parsed
  )

  pmin(pmax(norm, 0), 1)
}

fallback_rank_score <- function(position) {
  c(0.75, 0.5, 0.25)[position]
}

wilson_ci_scalar <- function(successes, total, conf = 0.95) {
  if (is.na(total) || total <= 0 || is.na(successes)) {
    return(c(NA_real_, NA_real_))
  }
  z <- qnorm(1 - (1 - conf) / 2)
  phat <- successes / total
  denom <- 1 + z^2 / total
  center <- (phat + z^2 / (2 * total)) / denom
  half_width <- (z / denom) * sqrt((phat * (1 - phat) / total) + z^2 / (4 * total^2))
  c(max(0, center - half_width), min(1, center + half_width))
}

compute_confusion_metrics <- function(df) {
  tp <- sum(df$truth_positive & df$pred_positive, na.rm = TRUE)
  fn <- sum(df$truth_positive & !df$pred_positive, na.rm = TRUE)
  tn <- sum(!df$truth_positive & !df$pred_positive, na.rm = TRUE)
  fp <- sum(!df$truth_positive & df$pred_positive, na.rm = TRUE)

  sens_denom <- tp + fn
  spec_denom <- tn + fp
  ppv_denom <- tp + fp
  npv_denom <- tn + fn

  sensitivity <- if (sens_denom > 0) tp / sens_denom else NA_real_
  specificity <- if (spec_denom > 0) tn / spec_denom else NA_real_
  ppv <- if (ppv_denom > 0) tp / ppv_denom else NA_real_
  npv <- if (npv_denom > 0) tn / npv_denom else NA_real_
  accuracy <- (tp + tn) / nrow(df)

  sens_ci <- wilson_ci_scalar(tp, sens_denom)
  spec_ci <- wilson_ci_scalar(tn, spec_denom)
  ppv_ci <- wilson_ci_scalar(tp, ppv_denom)
  npv_ci <- wilson_ci_scalar(tn, npv_denom)
  acc_ci <- wilson_ci_scalar(tp + tn, nrow(df))

  tibble(
    n = nrow(df),
    positive_n = sens_denom,
    negative_n = spec_denom,
    tp = tp,
    fn = fn,
    tn = tn,
    fp = fp,
    sensitivity = sensitivity,
    sensitivity_ci_low = sens_ci[1],
    sensitivity_ci_high = sens_ci[2],
    specificity = specificity,
    specificity_ci_low = spec_ci[1],
    specificity_ci_high = spec_ci[2],
    ppv = ppv,
    ppv_ci_low = ppv_ci[1],
    ppv_ci_high = ppv_ci[2],
    npv = npv,
    npv_ci_low = npv_ci[1],
    npv_ci_high = npv_ci[2],
    accuracy = accuracy,
    accuracy_ci_low = acc_ci[1],
    accuracy_ci_high = acc_ci[2]
  )
}

roc_pr_curve_metrics <- function(truth_positive, score) {
  df <- tibble(truth_positive = as.integer(truth_positive), score = score) %>%
    filter(!is.na(score))

  n_pos <- sum(df$truth_positive == 1)
  n_neg <- sum(df$truth_positive == 0)

  if (nrow(df) == 0 || n_pos == 0 || n_neg == 0) {
    return(list(
      auroc = NA_real_,
      auprc = NA_real_,
      roc_curve = tibble(fpr = NA_real_, tpr = NA_real_),
      pr_curve = tibble(recall = NA_real_, precision = NA_real_)
    ))
  }

  ranks <- rank(df$score, ties.method = "average")
  auroc <- (sum(ranks[df$truth_positive == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)

  thresholds <- sort(unique(df$score), decreasing = TRUE)
  roc_points <- map_dfr(thresholds, function(thr) {
    pred <- df$score >= thr
    tp <- sum(pred & df$truth_positive == 1)
    fp <- sum(pred & df$truth_positive == 0)
    tibble(
      threshold = thr,
      tpr = tp / n_pos,
      fpr = fp / n_neg
    )
  }) %>%
    arrange(fpr, tpr) %>%
    bind_rows(tibble(threshold = Inf, tpr = 0, fpr = 0), ., tibble(threshold = -Inf, tpr = 1, fpr = 1)) %>%
    distinct(fpr, tpr, .keep_all = TRUE) %>%
    arrange(fpr, tpr)

  pr_points <- map_dfr(thresholds, function(thr) {
    pred <- df$score >= thr
    tp <- sum(pred & df$truth_positive == 1)
    fp <- sum(pred & df$truth_positive == 0)
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 1
    recall <- tp / n_pos
    tibble(threshold = thr, recall = recall, precision = precision)
  }) %>%
    arrange(recall, desc(precision)) %>%
    distinct(recall, .keep_all = TRUE)

  pr_curve <- bind_rows(
    tibble(threshold = Inf, recall = 0, precision = 1),
    pr_points
  ) %>%
    arrange(recall)

  recall_diff <- diff(pr_curve$recall)
  auprc <- sum(pr_curve$precision[-1] * recall_diff)

  list(
    auroc = auroc,
    auprc = auprc,
    roc_curve = roc_points %>% select(fpr, tpr),
    pr_curve = pr_curve %>% select(recall, precision)
  )
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
  tibble(statistic = z_value, p_value = 2 * pnorm(abs(z_value), lower.tail = FALSE))
}

mcnemar_exact <- function(df, app_1, app_2) {
  b <- sum(df$app == app_1 & FALSE) # placeholder to keep linter quiet
}

pairwise_mcnemar_exact <- function(wide_df, app_1, app_2) {
  x <- wide_df[[app_1]]
  y <- wide_df[[app_2]]
  discordant_10 <- sum(x == 1 & y == 0, na.rm = TRUE)
  discordant_01 <- sum(x == 0 & y == 1, na.rm = TRUE)
  discordant_total <- discordant_10 + discordant_01
  p_value <- if (discordant_total > 0) {
    binom.test(min(discordant_10, discordant_01), discordant_total, p = 0.5)$p.value
  } else {
    NA_real_
  }
  tibble(
    app_1 = app_1,
    app_2 = app_2,
    discordant_10 = discordant_10,
    discordant_01 = discordant_01,
    p_value = p_value
  )
}

cochran_q_test <- function(wide_df) {
  mat <- as.matrix(wide_df[, unname(apps)])
  n <- nrow(mat)
  k <- ncol(mat)
  col_sums <- colSums(mat)
  row_sums <- rowSums(mat)
  total_success <- sum(col_sums)
  numerator <- (k - 1) * (k * sum(col_sums^2) - total_success^2)
  denominator <- k * total_success - sum(row_sums^2)
  if (denominator == 0) {
    return(tibble(statistic = NA_real_, df = k - 1, p_value = NA_real_))
  }
  q_value <- numerator / denominator
  tibble(
    statistic = q_value,
    df = k - 1,
    p_value = pchisq(q_value, df = k - 1, lower.tail = FALSE)
  )
}

message("Reading results from: ", results_path)

truth_diag <- read_csv(truth_diag_path, show_col_types = FALSE) %>%
  mutate(grouped_target = group_truth_disease(disease_category))

grouped_diag_terms <- truth_diag %>%
  group_by(grouped_target) %>%
  summarise(
    positive_terms = list(sort(unique(unlist(map(correct_diagnosis_classifications, split_terms))))),
    .groups = "drop"
  )

diag_term_lookup <- setNames(grouped_diag_terms$positive_terms, grouped_diag_terms$grouped_target)

truth_benmal <- read_csv(truth_benmal_path, show_col_types = FALSE) %>%
  mutate(phrases_list = map(phrases, split_terms))

benmal_lookup <- truth_benmal %>%
  select(benign_or_malignant, phrases_list) %>%
  tidyr::unnest_longer(phrases_list, values_to = "diagnosis") %>%
  transmute(diagnosis, benign_or_malignant) %>%
  distinct(diagnosis, .keep_all = TRUE)

raw <- read_csv(results_path, show_col_types = FALSE) %>%
  filter(!is.na(label), str_squish(label) != "") %>%
  mutate(
    fitzpatrick_scale = as.character(fitzpatrick_scale),
    fitzpatrick_group_3 = fitz_group_3(fitzpatrick_scale),
    truth_label_grouped = group_truth_disease(str_squish(label)),
    truth_malignancy = normalize_benmal(three_partition_label)
  ) %>%
  filter(truth_label_grouped %in% c("benign nevus", "melanoma spectrum", "basal cell carcinoma", "squamous cell carcinoma"))

imputation_audit <- imap_dfr(apps, function(app_name, app_code) {
  dx_cols <- paste0(app_code, "_dx_", 1:3)
  score_cols <- paste0(app_code, "_%_", 1:3)

  map_dfr(seq_len(nrow(raw)), function(i) {
    diagnoses <- map_chr(dx_cols, ~ str_squish(as.character(raw[[.x]][i])))
    scores <- map_dbl(score_cols, ~ parse_score(app_code, raw[[.x]][i]))
    target_terms <- list(
      malignant_vs_benign = benmal_lookup$diagnosis[
        benmal_lookup$benign_or_malignant == "malignant"
      ],
      melanoma_vs_rest = diag_term_lookup[["melanoma spectrum"]],
      basal_cell_carcinoma_vs_rest = diag_term_lookup[["basal cell carcinoma"]],
      squamous_cell_carcinoma_vs_rest = diag_term_lookup[["squamous cell carcinoma"]]
    )

    imap_dfr(target_terms, function(terms, target_name) {
      positive_idx <- which(diagnoses %in% terms)
      imputed_ranks <- positive_idx[is.na(scores[positive_idx])]
      tibble(
        target = target_name,
        imputed_differential_n = length(imputed_ranks),
        image_had_imputation = length(imputed_ranks) > 0
      )
    })
  }) %>%
    group_by(target) %>%
    summarise(
      imputed_differential_n = sum(imputed_differential_n),
      images_with_imputation_n = sum(image_had_imputation),
      .groups = "drop"
    ) %>%
    mutate(app = app_name, .before = 1)
})

write_csv(
  imputation_audit,
  file.path(output_dir, "confidence_score_imputation_counts.csv")
)

case_level <- imap_dfr(apps, function(app_name, app_code) {
  dx_cols <- paste0(app_code, "_dx_", 1:3)
  score_cols <- paste0(app_code, "_%_", 1:3)

  target_rows <- map_dfr(seq_len(nrow(raw)), function(i) {
    row <- raw[i, ]
    diagnoses <- map_chr(dx_cols, ~ str_squish(as.character(row[[.x]])))
    scores <- map2_dbl(score_cols, seq_along(score_cols), function(col, idx) {
      parsed <- parse_score(app_code, row[[col]])
      if (is.na(parsed)) NA_real_ else parsed
    })

    dx1 <- diagnoses[1]
    dx1_malignancy <- benmal_lookup$benign_or_malignant[match(dx1, benmal_lookup$diagnosis)]

    target_defs <- list(
      malignant_vs_benign = list(
        truth_positive = row$truth_malignancy == "malignant",
        pred_positive = !is.na(dx1_malignancy) && dx1_malignancy == "malignant",
        score = {
          positive_idx <- which(diagnoses %in% benmal_lookup$diagnosis[benmal_lookup$benign_or_malignant == "malignant"])
          if (length(positive_idx) == 0) {
            0
          } else {
            score_candidates <- map2_dbl(positive_idx, positive_idx, function(pos, idx) {
              if (!is.na(scores[pos])) scores[pos] else fallback_rank_score(pos)
            })
            max(score_candidates)
          }
        }
      ),
      melanoma_vs_rest = list(
        truth_positive = row$truth_label_grouped == "melanoma spectrum",
        pred_positive = dx1 %in% diag_term_lookup[["melanoma spectrum"]],
        score = {
          positive_idx <- which(diagnoses %in% diag_term_lookup[["melanoma spectrum"]])
          if (length(positive_idx) == 0) 0 else max(map2_dbl(positive_idx, positive_idx, function(pos, idx) if (!is.na(scores[pos])) scores[pos] else fallback_rank_score(pos)))
        }
      ),
      basal_cell_carcinoma_vs_rest = list(
        truth_positive = row$truth_label_grouped == "basal cell carcinoma",
        pred_positive = dx1 %in% diag_term_lookup[["basal cell carcinoma"]],
        score = {
          positive_idx <- which(diagnoses %in% diag_term_lookup[["basal cell carcinoma"]])
          if (length(positive_idx) == 0) 0 else max(map2_dbl(positive_idx, positive_idx, function(pos, idx) if (!is.na(scores[pos])) scores[pos] else fallback_rank_score(pos)))
        }
      ),
      squamous_cell_carcinoma_vs_rest = list(
        truth_positive = row$truth_label_grouped == "squamous cell carcinoma",
        pred_positive = dx1 %in% diag_term_lookup[["squamous cell carcinoma"]],
        score = {
          positive_idx <- which(diagnoses %in% diag_term_lookup[["squamous cell carcinoma"]])
          if (length(positive_idx) == 0) 0 else max(map2_dbl(positive_idx, positive_idx, function(pos, idx) if (!is.na(scores[pos])) scores[pos] else fallback_rank_score(pos)))
        }
      )
    )

    imap_dfr(target_defs, function(def, target_name) {
      tibble(
        md5hash = row$md5hash,
        fitzpatrick_scale = row$fitzpatrick_scale,
        fitzpatrick_group_3 = row$fitzpatrick_group_3,
        truth_label_grouped = row$truth_label_grouped,
        truth_malignancy = row$truth_malignancy,
        app = app_name,
        app_code = app_code,
        target = target_name,
        truth_positive = as.logical(def$truth_positive),
        pred_positive = as.logical(def$pred_positive),
        score = as.numeric(def$score),
        correct = as.integer(def$truth_positive == def$pred_positive)
      )
    })
  })
})

write_csv(case_level, file.path(output_dir, "case_level_binary_predictions_and_scores.csv"))

overall_metrics <- case_level %>%
  group_by(target, app) %>%
  group_modify(~ compute_confusion_metrics(.x)) %>%
  ungroup() %>%
  mutate(stratum = "Overall")

tone_metrics <- case_level %>%
  group_by(target, app, fitzpatrick_group_3) %>%
  group_modify(~ compute_confusion_metrics(.x)) %>%
  ungroup() %>%
  rename(stratum = fitzpatrick_group_3)

all_metrics <- bind_rows(overall_metrics, tone_metrics)
write_csv(all_metrics, file.path(output_dir, "binary_confusion_metrics_overall_and_by_tone.csv"))

discrimination <- case_level %>%
  group_by(target, app) %>%
  group_modify(~ {
    curves <- roc_pr_curve_metrics(.x$truth_positive, .x$score)
    tibble(
      n = nrow(.x),
      positive_n = sum(.x$truth_positive),
      negative_n = sum(!.x$truth_positive),
      auroc = curves$auroc,
      auprc = curves$auprc
    )
  }) %>%
  ungroup() %>%
  mutate(stratum = "Overall")

discrimination_tone <- case_level %>%
  group_by(target, app, fitzpatrick_group_3) %>%
  group_modify(~ {
    curves <- roc_pr_curve_metrics(.x$truth_positive, .x$score)
    tibble(
      n = nrow(.x),
      positive_n = sum(.x$truth_positive),
      negative_n = sum(!.x$truth_positive),
      auroc = curves$auroc,
      auprc = curves$auprc
    )
  }) %>%
  ungroup() %>%
  rename(stratum = fitzpatrick_group_3)

write_csv(bind_rows(discrimination, discrimination_tone), file.path(output_dir, "discrimination_metrics_overall_and_by_tone.csv"))

roc_curve_data <- case_level %>%
  filter(target == "malignant_vs_benign") %>%
  group_by(app) %>%
  group_modify(~ {
    curves <- roc_pr_curve_metrics(.x$truth_positive, .x$score)
    curves$roc_curve
  }) %>%
  ungroup()

pr_curve_data <- case_level %>%
  filter(target == "malignant_vs_benign") %>%
  group_by(app) %>%
  group_modify(~ {
    curves <- roc_pr_curve_metrics(.x$truth_positive, .x$score)
    curves$pr_curve
  }) %>%
  ungroup()

write_csv(roc_curve_data, file.path(output_dir, "malignant_vs_benign_roc_curve_points.csv"))
write_csv(pr_curve_data, file.path(output_dir, "malignant_vs_benign_pr_curve_points.csv"))

between_app_pairwise <- case_level %>%
  group_by(target) %>%
  group_modify(~ {
    wide_df <- .x %>%
      select(md5hash, app, correct) %>%
      distinct() %>%
      pivot_wider(names_from = app, values_from = correct)

    pairwise <- map_dfr(combn(unname(apps), 2, simplify = FALSE), function(pair) {
      pairwise_mcnemar_exact(wide_df, pair[1], pair[2])
    }) %>%
      mutate(p_value_holm = p.adjust(p_value, method = "holm"))

    q_test <- cochran_q_test(wide_df) %>%
      mutate(app_1 = NA_character_, app_2 = NA_character_, discordant_10 = NA_integer_, discordant_01 = NA_integer_, p_value_holm = NA_real_, test = "Cochran Q")

    pairwise %>% mutate(test = "McNemar exact") %>%
      bind_rows(q_test) %>%
      select(test, everything())
  }) %>%
  ungroup()

write_csv(between_app_pairwise, file.path(output_dir, "between_app_tests_overall.csv"))

tone_correctness_tests <- case_level %>%
  group_by(target, app) %>%
  group_modify(~ {
    tab <- table(.x$fitzpatrick_group_3, .x$correct)
    test_type <- safe_test_name(tab)
    test_obj <- if (test_type == "fisher") fisher.test(tab) else chisq.test(tab)

    counts <- .x %>%
      group_by(fitzpatrick_group_3) %>%
      summarise(correct_n = sum(correct), n = n(), .groups = "drop") %>%
      complete(fitzpatrick_group_3 = c("I-II", "III-IV", "V-VI"), fill = list(correct_n = 0, n = 0)) %>%
      arrange(fitzpatrick_group_3)
    trend <- cochran_armitage_test(counts$correct_n, counts$n, 1:3)

    tibble(
      test = if_else(test_type == "fisher", "Fisher exact", "Chi-square"),
      statistic = if (is.null(test_obj$statistic)) NA_real_ else unname(test_obj$statistic),
      df = if (is.null(test_obj$parameter)) NA_real_ else unname(test_obj$parameter),
      p_value = test_obj$p.value,
      ordinal_z = trend$statistic,
      ordinal_p_value = trend$p_value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_value_holm = p.adjust(p_value, method = "holm"),
    ordinal_p_value_holm = p.adjust(ordinal_p_value, method = "holm")
  )

write_csv(tone_correctness_tests, file.path(output_dir, "within_app_tone_tests_on_overall_correctness.csv"))

sensitivity_tone_tests <- case_level %>%
  filter(truth_positive) %>%
  group_by(target, app) %>%
  group_modify(~ {
    tab <- table(.x$fitzpatrick_group_3, .x$pred_positive)
    test_type <- safe_test_name(tab)
    test_obj <- if (test_type == "fisher") fisher.test(tab) else chisq.test(tab)
    counts <- .x %>%
      group_by(fitzpatrick_group_3) %>%
      summarise(success_n = sum(pred_positive), n = n(), .groups = "drop") %>%
      complete(fitzpatrick_group_3 = c("I-II", "III-IV", "V-VI"), fill = list(success_n = 0, n = 0)) %>%
      arrange(fitzpatrick_group_3)
    trend <- cochran_armitage_test(counts$success_n, counts$n, 1:3)
    tibble(
      test = if_else(test_type == "fisher", "Fisher exact", "Chi-square"),
      statistic = if (is.null(test_obj$statistic)) NA_real_ else unname(test_obj$statistic),
      df = if (is.null(test_obj$parameter)) NA_real_ else unname(test_obj$parameter),
      p_value = test_obj$p.value,
      ordinal_z = trend$statistic,
      ordinal_p_value = trend$p_value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_value_holm = p.adjust(p_value, method = "holm"),
    ordinal_p_value_holm = p.adjust(ordinal_p_value, method = "holm")
  )

write_csv(sensitivity_tone_tests, file.path(output_dir, "within_app_tone_tests_on_sensitivity.csv"))

specificity_tone_tests <- case_level %>%
  filter(!truth_positive) %>%
  mutate(pred_negative = !pred_positive) %>%
  group_by(target, app) %>%
  group_modify(~ {
    tab <- table(.x$fitzpatrick_group_3, .x$pred_negative)
    test_type <- safe_test_name(tab)
    test_obj <- if (test_type == "fisher") fisher.test(tab) else chisq.test(tab)
    counts <- .x %>%
      group_by(fitzpatrick_group_3) %>%
      summarise(success_n = sum(pred_negative), n = n(), .groups = "drop") %>%
      complete(fitzpatrick_group_3 = c("I-II", "III-IV", "V-VI"), fill = list(success_n = 0, n = 0)) %>%
      arrange(fitzpatrick_group_3)
    trend <- cochran_armitage_test(counts$success_n, counts$n, 1:3)
    tibble(
      test = if_else(test_type == "fisher", "Fisher exact", "Chi-square"),
      statistic = if (is.null(test_obj$statistic)) NA_real_ else unname(test_obj$statistic),
      df = if (is.null(test_obj$parameter)) NA_real_ else unname(test_obj$parameter),
      p_value = test_obj$p.value,
      ordinal_z = trend$statistic,
      ordinal_p_value = trend$p_value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_value_holm = p.adjust(p_value, method = "holm"),
    ordinal_p_value_holm = p.adjust(ordinal_p_value, method = "holm")
  )

write_csv(specificity_tone_tests, file.path(output_dir, "within_app_tone_tests_on_specificity.csv"))

disease_sensitivity_comparison <- case_level %>%
  filter(target %in% c("melanoma_vs_rest", "basal_cell_carcinoma_vs_rest", "squamous_cell_carcinoma_vs_rest"), truth_positive) %>%
  group_by(app) %>%
  group_modify(~ {
    tab <- table(.x$target, .x$pred_positive)
    test_type <- safe_test_name(tab)
    test_obj <- if (test_type == "fisher") fisher.test(tab) else chisq.test(tab)
    tibble(
      test = if_else(test_type == "fisher", "Fisher exact", "Chi-square"),
      statistic = if (is.null(test_obj$statistic)) NA_real_ else unname(test_obj$statistic),
      df = if (is.null(test_obj$parameter)) NA_real_ else unname(test_obj$parameter),
      p_value = test_obj$p.value
    )
  }) %>%
  ungroup()

write_csv(disease_sensitivity_comparison, file.path(output_dir, "within_app_disease_sensitivity_comparison_tests.csv"))

# Plots ---------------------------------------------------------------------
metrics_with_ci_long <- bind_rows(
  all_metrics %>%
    transmute(target, app, stratum, metric = "sensitivity", value = sensitivity, low = sensitivity_ci_low, high = sensitivity_ci_high),
  all_metrics %>%
    transmute(target, app, stratum, metric = "specificity", value = specificity, low = specificity_ci_low, high = specificity_ci_high),
  all_metrics %>%
    transmute(target, app, stratum, metric = "ppv", value = ppv, low = ppv_ci_low, high = ppv_ci_high),
  all_metrics %>%
    transmute(target, app, stratum, metric = "npv", value = npv, low = npv_ci_low, high = npv_ci_high)
) %>%
  mutate(
    app = factor(app, levels = app_levels),
    target = factor(target, levels = targets, labels = target_labels[targets]),
    metric = factor(metric, levels = c("sensitivity", "specificity", "ppv", "npv"),
                    labels = metric_labels[c("sensitivity", "specificity", "ppv", "npv")]),
    stratum = factor(stratum, levels = c("Overall", "I-II", "III-IV", "V-VI"))
  )

plot_overall_metrics <- ggplot(
  metrics_with_ci_long %>% filter(target == "Malignant vs Benign", stratum == "Overall"),
  aes(x = metric, y = value, color = app, group = app)
) +
  geom_point(position = position_dodge(width = 0.45), size = 2.8) +
  geom_errorbar(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.45),
    width = 0.18,
    linewidth = 0.6
  ) +
  scale_color_manual(values = app_palette) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Malignant vs Benign Screening Metrics by App",
    x = NULL,
    y = "Estimate (95% CI)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(
  filename = file.path(output_dir, "malignant_vs_benign_overall_metrics.png"),
  plot = plot_overall_metrics,
  width = 9.5,
  height = 6,
  dpi = 300
)

plot_tone_metrics <- ggplot(
  metrics_with_ci_long %>% filter(target == "Malignant vs Benign", stratum != "Overall"),
  aes(x = stratum, y = value, color = app, group = app)
) +
  geom_point(position = position_dodge(width = 0.45), size = 2.4) +
  geom_errorbar(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.45),
    width = 0.18,
    linewidth = 0.55
  ) +
  facet_wrap(~ metric, ncol = 2) +
  scale_color_manual(values = app_palette) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Malignant vs Benign Metrics by Skin Tone Group",
    x = NULL,
    y = "Estimate (95% CI)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(
  filename = file.path(output_dir, "malignant_vs_benign_metrics_by_tone.png"),
  plot = plot_tone_metrics,
  width = 11,
  height = 8,
  dpi = 300
)

plot_roc <- ggplot(roc_curve_data, aes(x = fpr, y = tpr, color = app)) +
  geom_path(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
  coord_equal() +
  scale_color_manual(values = app_palette) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "ROC Curves for Malignant vs Benign Classification",
    x = "False-positive rate",
    y = "True-positive rate",
    color = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(output_dir, "malignant_vs_benign_roc_curves.png"),
  plot = plot_roc,
  width = 8,
  height = 6,
  dpi = 300
)

plot_pr <- ggplot(pr_curve_data, aes(x = recall, y = precision, color = app)) +
  geom_path(linewidth = 1) +
  scale_color_manual(values = app_palette) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Precision-Recall Curves for Malignant vs Benign Classification",
    x = "Recall",
    y = "Precision",
    color = NULL
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(output_dir, "malignant_vs_benign_pr_curves.png"),
  plot = plot_pr,
  width = 8,
  height = 6,
  dpi = 300
)

discrimination_heatmap <- discrimination %>%
  filter(target != "malignant_vs_benign") %>%
  select(target, app, auroc, auprc) %>%
  pivot_longer(cols = c(auroc, auprc), names_to = "metric", values_to = "value") %>%
  mutate(
    target = factor(target, levels = targets[targets != "malignant_vs_benign"],
                    labels = target_labels[targets[targets != "malignant_vs_benign"]]),
    app = factor(app, levels = app_levels),
    metric = factor(metric, levels = c("auroc", "auprc"), labels = metric_labels[c("auroc", "auprc")])
  )

plot_discrimination_heatmap <- ggplot(discrimination_heatmap, aes(x = app, y = target, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  facet_wrap(~ metric) +
  scale_fill_gradient(low = "#f2e8cf", high = "#386641", na.value = "grey90") +
  labs(
    title = "Disease-vs-Rest Discrimination Within App",
    x = NULL,
    y = NULL,
    fill = "Value"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_dir, "disease_vs_rest_discrimination_heatmap.png"),
  plot = plot_discrimination_heatmap,
  width = 12,
  height = 5.8,
  dpi = 300
)

disease_vs_rest_metric_plot_data <- metrics_with_ci_long %>%
  filter(target != "Malignant vs Benign", stratum == "Overall")

plot_disease_vs_rest_metrics <- ggplot(
  disease_vs_rest_metric_plot_data,
  aes(x = metric, y = value, color = app, group = app)
) +
  geom_point(position = position_dodge(width = 0.45), size = 2.4) +
  geom_errorbar(
    aes(ymin = low, ymax = high),
    position = position_dodge(width = 0.45),
    width = 0.16,
    linewidth = 0.55
  ) +
  facet_wrap(~ target, ncol = 1) +
  scale_color_manual(values = app_palette) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Disease-vs-Rest Screening Metrics by App",
    x = NULL,
    y = "Estimate (95% CI)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_dir, "disease_vs_rest_screening_metrics.png"),
  plot = plot_disease_vs_rest_metrics,
  width = 10.5,
  height = 11,
  dpi = 300
)

message("Binary screening analysis complete.")
message("Outputs written to: ", output_dir)
