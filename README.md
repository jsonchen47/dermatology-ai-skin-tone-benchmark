# Consumer-Facing Dermatology AI Benchmark

This repository is the reproducibility package for **Diagnostic Performance of Consumer-Facing Artificial Intelligence Applications for Skin Cancer Across Fitzpatrick Skin Tones**.

The primary cohort contains 231 Fitzpatrick17k images evaluated by AI Skin Scanner, Model Dermatol, ChatGPT, and Claude. The repository includes cohort identifiers, source-derived labels, per-image application outputs, adjudication dictionaries, analysis-ready outcomes, R scripts, result tables, and publication figures. Source images are not redistributed.

## Cohorts

| Cohort | Disease composition | Grouped skin-tone composition |
|---|---:|---:|
| Fitzpatrick17k primary cohort | BCC 58; benign nevus 53; melanoma spectrum 60; SCC 60 | I-II 81; III-IV 80; V-VI 70 |
| DDI external cohort | See `data/ddi/ddi_cohort_252.csv` | I-II 84; III-IV 84; V-VI 84 |

The primary cohort is balanced as closely as the eligible Fitzpatrick17k pool permitted across the four disease groups and three grouped Fitzpatrick strata. It is enriched for malignant lesions and is not intended to estimate screening-population predictive values.

## Repository Structure

- `data/fitzpatrick/`: definitive 231-image manifest, per-image application outputs, analysis-ready correctness indicators, and diagnosis dictionaries.
- `data/ddi/`: 252-image DDI manifest and ChatGPT/Claude classifications used for the external analysis.
- `analysis/`: ordered R scripts, input validation, dependency helper, and the complete rerun script.
- `results/`: regenerated tables and intermediate analysis outputs.
- `figures/`: publication-ready PNG and vector PDF figures.
- `docs/`: data dictionary, analysis decisions, model prompt record, release checklist, and manuscript Data Sharing Statement.
- `environment/`: R and package-version records from the verified run.

## Reproduce the Analysis

Run from any directory after cloning the repository:

```bash
Rscript analysis/00_validate_inputs.R
Rscript analysis/run_all.R
```

If required packages are missing:

```bash
Rscript analysis/install_packages.R
```

The verified run used R 4.5.3. `analysis/06_publication_figures.R` uses a fixed random seed and 2,000 class-stratified bootstrap samples for AUROC and AUPRC confidence intervals. A full rerun may take several minutes.

## Analysis Map

1. `00_validate_inputs.R` confirms sample counts, unique identifiers, tone distributions, and agreement among deposited analysis files.
2. `01_primary_accuracy.R` creates descriptive accuracy summaries and supporting analyses.
3. `02_binary_screening.R` creates malignant-vs-benign and disease-vs-rest confusion metrics, exact binary comparisons, and within-application AUROC/AUPRC inputs.
4. `03_clustered_logistic_regression.R` fits adjusted logistic models with image-clustered robust standard errors.
5. `04_endpoint_family_tone_analysis.R` performs the four-application, family-level, disease-level, and grouped-tone analyses used for the principal figures.
6. `05_broad_vs_specific_analysis.R` compares broad diagnostic scoring with exact disease-level top-1 scoring.
7. `06_publication_figures.R` generates the publication figures and bootstrap confidence intervals.
8. `07_ddi_external_analysis.R` reproduces the secondary DDI evaluation.
9. `08_power_context.R` reproduces the sample-count and minimum-detectable-difference context tables.

## Source Datasets

- Fitzpatrick17k: https://github.com/mattgroh/fitzpatrick17k
- Diverse Dermatology Images (DDI): https://aimi.stanford.edu/datasets/ddi-diverse-dermatology-images

Download and use source images only under the original dataset terms. The MD5-hash list in `data/fitzpatrick/cohort_manifest_231.csv` is the definitive record of primary-cohort membership.

## Citation and Archiving

Citation metadata are provided in `CITATION.cff`. Before manuscript submission, create a versioned GitHub release, archive that release with Zenodo, and replace the placeholders in `docs/DATA_SHARING_STATEMENT.md` with the public repository URL and version-specific DOI.

## License

Analysis code is licensed under the MIT License. See `DATA_USE_NOTICE.md` for the separate terms and source-dataset limitations that apply to tabular data and derived outputs.
