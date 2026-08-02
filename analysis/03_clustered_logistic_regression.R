#!/usr/bin/env Rscript

# Image-clustered logistic regression for the Fitzpatrick 17K analysis.
# Each image contributes one prediction from each of four applications.
# Coefficients are estimated by ordinary logistic regression, while inference
# uses a cluster-robust sandwich covariance matrix with md5hash as the cluster.
# This is equivalent to an independence-working-correlation GEE for coefficients
# and robust standard errors.

required_packages <- c("dplyr", "tidyr", "readr", "stringr", "forcats", "purrr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(purrr)
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
  unset = file.path(repo_dir, "results", "clustered_regression")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

apps <- c(
  a = "AI Skin Scanner",
  m = "Model Dermatol",
  c = "ChatGPT",
  cl = "Claude"
)

outcomes <- c(
  dx_1_correct = "Top-1 diagnosis accuracy",
  any_dx_correct = "Top-3 diagnosis accuracy",
  dx_1_malignancy_correct = "Benign/malignant accuracy"
)

raw <- read_csv(input_path, show_col_types = FALSE) %>%
  mutate(
    fitzpatrick_scale = factor(as.character(fitzpatrick_scale), levels = as.character(1:6)),
    label = str_squish(label),
    label = case_when(
      label %in% c("congenital nevus", "halo nevus") ~ "benign nevus",
      label %in% c("malignant melanoma", "melanoma") ~ "melanoma spectrum",
      TRUE ~ label
    ),
    label = factor(
      label,
      levels = c(
        "melanoma spectrum",
        "squamous cell carcinoma",
        "basal cell carcinoma",
        "benign nevus"
      )
    )
  )

long <- raw %>%
  pivot_longer(
    cols = matches("^(a|m|c|cl)_(dx_1_correct|any_dx_correct|dx_1_malignancy_correct)$"),
    names_to = c("app_code", "outcome_code"),
    names_pattern = "^(a|m|c|cl)_(dx_1_correct|any_dx_correct|dx_1_malignancy_correct)$",
    values_to = "correct"
  ) %>%
  mutate(
    app = factor(
      unname(apps[app_code]),
      levels = c("ChatGPT", "AI Skin Scanner", "Model Dermatol", "Claude")
    ),
    outcome = unname(outcomes[outcome_code]),
    correct = as.integer(correct)
  )

cluster_robust_fit <- function(data) {
  fit <- glm(
    correct ~ app + fitzpatrick_scale + label,
    family = binomial(),
    data = data
  )

  X <- model.matrix(fit)
  y <- model.response(model.frame(fit))
  mu <- fitted(fit)
  cluster <- data$md5hash

  # Bread for logistic maximum likelihood.
  w <- mu * (1 - mu)
  bread_inv <- solve(crossprod(X, X * as.numeric(w)))

  # Sum observation-level score contributions within image.
  score_obs <- X * as.numeric(y - mu)
  cluster_scores <- rowsum(score_obs, group = cluster, reorder = FALSE)
  meat <- crossprod(cluster_scores)

  n <- nrow(X)
  k <- ncol(X)
  g <- nrow(cluster_scores)
  correction <- (g / (g - 1)) * ((n - 1) / (n - k))
  vcov_cluster <- correction * bread_inv %*% meat %*% bread_inv

  beta <- coef(fit)
  robust_se <- sqrt(diag(vcov_cluster))
  ordinary_se <- sqrt(diag(vcov(fit)))

  coef_table <- tibble(
    term = names(beta),
    log_odds = unname(beta),
    odds_ratio = exp(log_odds),
    ordinary_se = unname(ordinary_se),
    ordinary_p = 2 * pnorm(abs(log_odds / ordinary_se), lower.tail = FALSE),
    cluster_robust_se = unname(robust_se),
    cluster_robust_z = log_odds / cluster_robust_se,
    cluster_robust_p = 2 * pnorm(abs(cluster_robust_z), lower.tail = FALSE),
    conf_low = exp(log_odds - qnorm(0.975) * cluster_robust_se),
    conf_high = exp(log_odds + qnorm(0.975) * cluster_robust_se),
    n_observations = n,
    n_images = g
  )

  groups <- list(
    application = grep("^app", names(beta)),
    fitzpatrick_skin_type = grep("^fitzpatrick_scale", names(beta)),
    disease = grep("^label", names(beta))
  )

  joint_tests <- imap_dfr(groups, function(idx, group_name) {
    b <- beta[idx]
    v <- vcov_cluster[idx, idx, drop = FALSE]
    statistic <- as.numeric(t(b) %*% solve(v, b))
    tibble(
      predictor = group_name,
      df = length(idx),
      chi_square = statistic,
      cluster_robust_p = pchisq(statistic, df = length(idx), lower.tail = FALSE),
      n_observations = n,
      n_images = g
    )
  })

  list(fit = fit, vcov = vcov_cluster, coefficients = coef_table, joint = joint_tests)
}

fits <- long %>%
  split(.$outcome) %>%
  map(cluster_robust_fit)

coefficient_results <- imap_dfr(fits, function(result, outcome_name) {
  result$coefficients %>% mutate(outcome = outcome_name, .before = 1)
}) %>%
  mutate(
    p_holm_within_endpoint = if_else(
      term == "(Intercept)",
      NA_real_,
      ave(
        if_else(term == "(Intercept)", NA_real_, cluster_robust_p),
        outcome,
        FUN = function(x) {
          result <- rep(NA_real_, length(x))
          keep <- !is.na(x)
          result[keep] <- p.adjust(x[keep], method = "holm")
          result
        }
      )
    )
  )

non_intercept <- coefficient_results$term != "(Intercept)"
coefficient_results$p_holm_all_endpoints <- NA_real_
coefficient_results$p_holm_all_endpoints[non_intercept] <- p.adjust(
  coefficient_results$cluster_robust_p[non_intercept],
  method = "holm"
)

joint_results <- imap_dfr(fits, function(result, outcome_name) {
  result$joint %>% mutate(outcome = outcome_name, .before = 1)
}) %>%
  group_by(outcome) %>%
  mutate(p_holm_within_endpoint = p.adjust(cluster_robust_p, method = "holm")) %>%
  ungroup() %>%
  mutate(p_holm_all_joint_tests = p.adjust(cluster_robust_p, method = "holm"))

model_summary <- imap_dfr(fits, function(result, outcome_name) {
  tibble(
    outcome = outcome_name,
    n_observations = nobs(result$fit),
    n_images = n_distinct(long$md5hash),
    applications_per_image = 4,
    covariance = "Image-clustered sandwich (CR1)",
    cluster_id = "md5hash"
  )
})

write_csv(
  coefficient_results,
  file.path(output_dir, "clustered_logistic_regression_coefficients.csv")
)
write_csv(
  joint_results,
  file.path(output_dir, "clustered_logistic_regression_joint_tests.csv")
)
write_csv(
  model_summary,
  file.path(output_dir, "clustered_logistic_regression_model_summary.csv")
)

capture.output(
  {
    cat("IMAGE-CLUSTERED LOGISTIC REGRESSION\n")
    cat("Input: data/fitzpatrick/fitzpatrick_analysis_231.csv\n")
    cat("Clusters: 231 images; four application predictions per image and endpoint.\n")
    cat("Reference levels: ChatGPT, Fitzpatrick 1, melanoma spectrum.\n\n")
    cat("JOINT CLUSTER-ROBUST WALD TESTS\n")
    print(joint_results, n = Inf)
    cat("\nCOEFFICIENTS WITH RAW CLUSTER-ROBUST P < .05\n")
    print(
      coefficient_results %>%
        filter(term != "(Intercept)", cluster_robust_p < 0.05) %>%
        select(
          outcome, term, odds_ratio, conf_low, conf_high,
          cluster_robust_p, p_holm_within_endpoint
        ),
      n = Inf
    )
  },
  file = file.path(output_dir, "clustered_logistic_regression_summary.txt")
)

message("Wrote clustered regression outputs to: ", output_dir)
