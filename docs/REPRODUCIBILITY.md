# Reproducibility Notes

## Definitive Primary Cohort

The definitive primary cohort is the 231-row file `data/fitzpatrick/cohort_manifest_231.csv`. It contains 58 BCC, 53 benign-nevus, 60 melanoma-spectrum, and 60 SCC images. The grouped Fitzpatrick counts are 81 I-II, 80 III-IV, and 70 V-VI. Lentigo maligna cases are not included in the deposited cohort.

The working file from which this package was assembled contained five malformed or non-evaluable rows in addition to the 231 valid cases. Those rows are not included in the deposited cohort or analysis files. The input-validation script confirms the row counts, unique identifiers, and exact identifier agreement across all primary analysis files.

## Reproduction Order

`analysis/run_all.R` validates inputs, runs the analyses in dependency order, regenerates tables, calculates bootstrap intervals, and copies final PNG/PDF figures to `figures/`. The broad-versus-specific script uses the deposited no-lentigo per-image outputs rather than the obsolete lentigo-containing working file.

## Confidence-Score Audit

`results/binary_screening/confidence_score_imputation_counts.csv` reports, by application and discrimination target, how many mapped differentials lacked a confidence value and therefore received the prespecified rank-based fallback score.

## Scope

The repository reproduces reported analyses from deposited application outputs. It does not rerun the four proprietary applications against source images and does not guarantee that future versions of those applications will reproduce the historical outputs.
