#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  normalizePath("analysis")
}
repo_dir <- normalizePath(file.path(script_dir, ".."))
old_working_dir <- setwd(repo_dir)
on.exit(setwd(old_working_dir), add = TRUE)

scripts <- c(
  "00_validate_inputs.R",
  "01_primary_accuracy.R",
  "02_binary_screening.R",
  "03_clustered_logistic_regression.R",
  "04_endpoint_family_tone_analysis.R",
  "05_broad_vs_specific_analysis.R",
  "07_ddi_external_analysis.R",
  "08_power_context.R",
  "06_publication_figures.R"
)

for (script in scripts) {
  message("\nRunning ", script, "...")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    file.path("analysis", script),
    stdout = "",
    stderr = ""
  )
  if (!identical(status, 0L)) {
    stop("Analysis stopped because ", script, " returned status ", status, ".")
  }
}

publication_dir <- file.path(repo_dir, "results", "publication")
main_dir <- file.path(repo_dir, "figures", "main")
supp_dir <- file.path(repo_dir, "figures", "supplementary")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

main_figures <- list.files(
  publication_dir,
  pattern = "^figure[1-5].*\\.(png|pdf)$",
  full.names = TRUE
)
supp_figures <- list.files(
  publication_dir,
  pattern = "^supplementary_figure.*\\.(png|pdf)$",
  full.names = TRUE
)
invisible(file.copy(main_figures, main_dir, overwrite = TRUE))
invisible(file.copy(supp_figures, supp_dir, overwrite = TRUE))

message("\nAll analyses completed. Publication figures are in: ", file.path(repo_dir, "figures"))
