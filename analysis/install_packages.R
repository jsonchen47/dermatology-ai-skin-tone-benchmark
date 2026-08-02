#!/usr/bin/env Rscript

required_packages <- c(
  "broom", "dplyr", "forcats", "ggplot2", "patchwork",
  "purrr", "readr", "scales", "stringr", "tidyr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) == 0) {
  message("All required packages are already installed.")
} else {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}
