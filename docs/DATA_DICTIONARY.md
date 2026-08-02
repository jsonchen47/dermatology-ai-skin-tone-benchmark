# Data Dictionary

## Application Prefixes

| Prefix | Application |
|---|---|
| `a` | AI Skin Scanner |
| `m` | Model Dermatol |
| `c` | ChatGPT |
| `cl` | Claude |

## Primary Cohort Files

### `cohort_manifest_231.csv`

| Field | Definition |
|---|---|
| `md5hash` | Fitzpatrick17k image identifier; definitive key for cohort membership. |
| `fitzpatrick_scale` | Scale AI Fitzpatrick skin-type label used in the analysis, from 1 to 6. |
| `grouped_fitzpatrick_tone` | Analysis group: I-II, III-IV, or V-VI. |
| `source_label` | Original disease label retained from Fitzpatrick17k. |
| `disease_group` | Four-group study label: BCC, benign nevus, melanoma spectrum, or SCC. |
| `malignancy` | Source three-partition designation used as benign/malignant truth. |

### `per_image_application_outputs_231.csv`

The file includes the cohort metadata above plus application outputs. For each application prefix:

| Suffix | Definition |
|---|---|
| `_score` | Application-reported overall dangerousness or malignancy score in its native scale; not used for cross-application AUROC comparison. |
| `_dx_1`, `_dx_2`, `_dx_3` | Ranked first, second, and third diagnostic differentials. |
| `_%_1`, `_%_2`, `_%_3` | Reported confidence associated with each ranked differential. Scales are normalized within the analysis script. |
| `_notes` | Application-generated notes retained during data collection. |

`fitzpatrick_centaur`, `nine_partition_label`, and `three_partition_label` are source-derived metadata fields. URLs and source images are intentionally omitted from this public analysis file.

### Analysis-Ready Files

`fitzpatrick_analysis_231.csv` contains binary correctness indicators for broad diagnostic scoring. `fitzpatrick_analysis_exact_top1_231.csv` contains the stricter exact top-1 indicators. The indicator suffixes are:

| Suffix | Definition |
|---|---|
| `_dx_1_correct` | Top-ranked diagnosis met the applicable phrase-dictionary rule. |
| `_any_dx_correct` | At least one of the top 3 diagnoses met the applicable phrase-dictionary rule. |
| `_dx_1_malignancy_correct` | Top-ranked diagnosis mapped to the correct benign/malignant class. |

### Dictionaries

`truth_correct_diagnosis.csv` defines accepted broad diagnostic phrases. `truth_correct_diagnosis_exact_top1.csv` defines the disease-level top-1 sensitivity analysis and records terms considered too broad. `truth_malignant_or_benign.csv` maps diagnostic phrases to binary malignancy status.

## DDI Files

`ddi_cohort_252.csv` contains the DDI identifier, filename, grouped skin tone (`12`, `34`, or `56`), source malignancy label, and disease. `chatgpt_outputs_252.csv` and `claude_outputs_252.csv` contain each model's diagnosis and benign/malignant classification for the same identifiers.

## Missing Values and Confidence Scores

Blank CSV cells represent unavailable or unreported values. For discrimination analyses, confidence values were normalized to 0 to 1. The score for a target was the highest confidence among top-3 differentials mapped to that target; when a mapped differential lacked confidence, rank-based values of 0.75, 0.50, and 0.25 were used for ranks 1, 2, and 3. Images with no mapped positive differential received a score of 0. AUROC/AUPRC are interpreted within applications because confidence scales were not calibrated across applications.

For the main malignant-vs-benign analysis, 3 Model Dermatol differentials across 3 images required rank-based confidence imputation; no differentials required imputation for AI Skin Scanner, ChatGPT, or Claude. The complete target-specific audit is in `results/binary_screening/confidence_score_imputation_counts.csv`.
