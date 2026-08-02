#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  normalizePath("analysis")
}
repo_dir <- normalizePath(file.path(script_dir, ".."))

fitz_dir <- file.path(repo_dir, "data", "fitzpatrick")
ddi_dir <- file.path(repo_dir, "data", "ddi")
output_path <- file.path(repo_dir, "results", "input_validation.txt")

outputs <- read.csv(
  file.path(fitz_dir, "per_image_application_outputs_231.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
manifest <- read.csv(
  file.path(fitz_dir, "cohort_manifest_231.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
analysis <- read.csv(
  file.path(fitz_dir, "fitzpatrick_analysis_231.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
exact <- read.csv(
  file.path(fitz_dir, "fitzpatrick_analysis_exact_top1_231.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
ddi <- read.csv(
  file.path(ddi_dir, "ddi_cohort_252.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(outputs) == 231L,
  nrow(manifest) == 231L,
  nrow(analysis) == 231L,
  nrow(exact) == 231L,
  length(unique(outputs$md5hash)) == 231L,
  setequal(outputs$md5hash, manifest$md5hash),
  setequal(outputs$md5hash, analysis$md5hash),
  setequal(outputs$md5hash, exact$md5hash),
  nrow(ddi) == 252L,
  length(unique(ddi$DDI_ID)) == 252L
)

expected_disease <- c(
  "basal cell carcinoma" = 58L,
  "benign nevus" = 53L,
  "melanoma spectrum" = 60L,
  "squamous cell carcinoma" = 60L
)
expected_tone <- c("I-II" = 81L, "III-IV" = 80L, "V-VI" = 70L)
expected_ddi_tone <- c("12" = 84L, "34" = 84L, "56" = 84L)

count_matches <- function(values, expected) {
  observed <- table(values)
  all(names(expected) %in% names(observed)) &&
    all(as.integer(observed[names(expected)]) == unname(expected)) &&
    sum(observed) == sum(expected)
}

stopifnot(
  count_matches(manifest$disease_group, expected_disease),
  count_matches(manifest$grouped_fitzpatrick_tone, expected_tone),
  count_matches(as.character(ddi$skin_tone), expected_ddi_tone)
)

data_files <- list.files(
  file.path(repo_dir, "data"),
  recursive = TRUE,
  full.names = TRUE
)
data_checksums <- tools::md5sum(data_files)
checksum_names <- sub(paste0("^", repo_dir, "/"), "", names(data_checksums))

lines <- c(
  "Input validation passed.",
  "",
  paste0("Primary rows: ", nrow(outputs)),
  paste0("Unique primary MD5 identifiers: ", length(unique(outputs$md5hash))),
  paste0("Primary disease counts: ", paste(names(table(manifest$disease_group)), table(manifest$disease_group), collapse = "; ")),
  paste0("Primary grouped-tone counts: ", paste(names(table(manifest$grouped_fitzpatrick_tone)), table(manifest$grouped_fitzpatrick_tone), collapse = "; ")),
  paste0("DDI rows: ", nrow(ddi)),
  paste0("DDI grouped-tone counts: ", paste(names(table(ddi$skin_tone)), table(ddi$skin_tone), collapse = "; ")),
  "",
  "Input MD5 checksums:",
  paste(checksum_names, unname(data_checksums), sep = "  ")
)

writeLines(lines, output_path)
cat(paste(lines[1:8], collapse = "\n"), "\n")
