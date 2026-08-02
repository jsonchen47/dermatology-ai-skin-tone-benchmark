# Analysis Decisions

## Paired Design

Every primary-cohort image was evaluated by all four applications. Cochran Q tests assess overall application differences for binary correctness endpoints. Pairwise application comparisons for the main diagnostic-accuracy endpoints use continuity-corrected McNemar tests; exact binomial McNemar tests are used in the binary screening script. Holm adjustments are applied within the prespecified comparison families implemented in each script.

## Application-Family Comparison

For each image and endpoint, the two smartphone applications receive a family score equal to their mean correctness (0, 0.5, or 1), and ChatGPT and Claude receive the analogous LLM-family score. A paired Wilcoxon signed-rank test compares the two image-level family scores. The three endpoint P values are Holm-adjusted together.

## Skin-Tone Analysis

Grouped associations use I-II, III-IV, and V-VI. Overall contingency-table tests assess any difference among the three strata, while ordinal trend tests assess a monotonic pattern. Holm correction is applied separately across the four applications for overall and ordinal primary tone tests.

## Adjusted Regression

The adjusted models include application, six-level Fitzpatrick skin type, and four-group disease category. Because each image contributes one observation per application, `03_clustered_logistic_regression.R` uses image-clustered robust standard errors. These models are exploratory and support endpoint-specific interpretation.

## Binary Screening and Discrimination

Sensitivity, specificity, PPV, NPV, and accuracy are calculated from malignant-vs-benign and disease-vs-rest confusion matrices. Wilson 95% confidence intervals are used for binary proportions. AUROC and AUPRC use the highest confidence among top-3 differentials mapped to the positive target, not the application's separate overall dangerousness score. Confidence scales differ by application, so AUROC/AUPRC are reported within applications and are not formally compared across applications.

The publication figure script uses 2,000 class-stratified bootstrap replicates and seed `20260725` for AUROC/AUPRC 95% confidence intervals.

## DDI Analysis

DDI is an external evaluation of ChatGPT and Claude, not a causal or direct estimate of the effect of dataset balancing. Cohorts differ in disease composition, image characteristics, and diagnostic scoring rules.
