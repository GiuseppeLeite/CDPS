# Integrated Profiling of Non-Apoptotic Cell Death Identifies SLC3A2-Mediated Ferroptosis Suppression as a Prognostic Determinant in Melanoma

This repository provides the final trained Cell Death Programmed Score (CDPS) model and an example script showing how to apply the model to an external melanoma cohort.

The CDPS was developed as described in the manuscript using genes associated with regulated non-apoptotic cell death pathways, including necroptosis, pyroptosis, and ferroptosis. Briefly, candidate genes were derived from curated cell death pathways, expression matrices were harmonized using Z-score normalization, and machine-learning survival models were trained and validated across independent melanoma cohorts.

Model development was performed using the Mime R package, a machine-learning framework for predictive model construction, visualization, and feature selection:

https://github.com/l-magnificence/Mime

This repository provides the final trained model object and a reproducible example showing how to calculate CDPS in an external dataset.

## Repository contents

```text
.
├── CDPS_code.R
├── model_rsf_ridge.rds
├── GSE91061_mime.rds
└── README.md
```

## Files

### `CDPS_code.R`

R script used to apply the trained CDPS model to the example dataset and generate an illustrative paired pre-treatment versus on-treatment plot for GSE91061.

### `model_rsf_ridge.rds`

Final trained CDPS model object.

This object was generated during the original model-development workflow using the Mime R package. It contains the final selected ridge-based model used for CDPS calculation.

The model object contains the fitted ridge model and the selected lambda value used for prediction:

```r
model_rsf_ridge$glmnet.fit
model_rsf_ridge$lambda.min
```

The genes used by the final model are extracted directly from the trained model object:

```r
rownames(model_rsf_ridge$glmnet.fit$beta)
```

### `GSE91061_mime.rds`

Example input dataset used to illustrate CDPS calculation and plotting.

This file is included only as an example of the expected structure. Users may apply the CDPS model to their own datasets, provided that the required model genes are present as columns.

## Input data format

The input dataset must be readable with `readRDS()` and convertible to a `data.frame`.

Rows must represent samples. Columns must include the genes required by the trained CDPS model.

Minimal structure required for CDPS calculation:

```text
Sample_ID   Gene1   Gene2   Gene3   Gene4   Gene5   Gene6   Gene7
S001        ...     ...     ...     ...     ...     ...     ...
S002        ...     ...     ...     ...     ...     ...     ...
S003        ...     ...     ...     ...     ...     ...     ...
```

The exact required gene symbols are not manually specified in the script. They are extracted from:

```r
rownames(model_rsf_ridge$glmnet.fit$beta)
```

All model genes must be present as columns in the input dataset.

The script stops if:

```text
- any model gene is missing;
- gene expression values are missing or non-numeric;
- any model gene has zero variance after Z-score scaling.
```

## Expression data requirements

Input expression values should already be normalized appropriately for the platform used.

For RNA-seq datasets, TPM or equivalent normalized expression values are recommended. For microarray datasets, platform-normalized expression values mapped to gene symbols can be used.

The script applies gene-wise Z-score normalization within the input dataset before calculating the CDPS, consistent with the harmonization strategy described in the Methods section of the manuscript.

The input data should therefore not be manually Z-scored before running the script unless the user intentionally modifies the script accordingly.

## Clinical or survival metadata

Clinical metadata are not required to calculate the CDPS.

Users may include any additional clinical or survival variables relevant to their own analysis, for example:

```text
OS_time
OS_status
DSS_time
DSS_status
DFS_time
DFS_status
PFS_time
PFS_status
Stage
Treatment
Response
```

These variables are not used by the CDPS prediction step itself. They can be used downstream to test associations between CDPS and survival, treatment response, clinical phenotype, or other outcomes.

## GSE91061 example metadata

The included `GSE91061_mime.rds` example contains additional columns used only for the illustrative Pre versus On-treatment plot.

For this example, the plotting section of `CDPS_code.R` expects:

```text
Patient
Group
Response
```

Expected coding in the example:

```text
Group:
  Pre = pre-treatment sample
  On  = on-treatment sample

Response:
  PRCR = responder
  other values = non-responder
```

These columns are specific to the GSE91061 example and are not required for CDPS calculation in other datasets.

## CDPS calculation

Expression values from the model genes are extracted in the exact order used by the trained ridge model. The expression matrix is then Z-score standardized by gene within the input dataset.

The CDPS is calculated using:

```r
predict(
  model_rsf_ridge$glmnet.fit,
  newx = expr_z,
  s = model_rsf_ridge$lambda.min,
  type = "link"
)
```

No gene imputation is performed.

## Required R packages

The script requires:

```r
install.packages(c("tidyverse", "glmnet"))
```

## Running the analysis

Run the script in R from the repository directory:

```r
source("CDPS_code.R")
```

## Output files

The script generates:

```text
GSE91061_mime_CDPS.rds
GSE91061_PRE_ON_ByResponse_CDPS.pdf
```

### `GSE91061_mime_CDPS.rds`

The example input dataset with an additional column:

```text
CDPS
```

### `GSE91061_PRE_ON_ByResponse_CDPS.pdf`

Illustrative paired pre-treatment versus on-treatment CDPS plot stratified by response group in GSE91061.

## Notes

This repository is intended to document and reproduce the application of the final trained CDPS model to an external melanoma cohort.

The full model-development strategy, including candidate gene selection, training cohorts, validation cohorts, feature-selection procedure, algorithm comparison, and survival-performance evaluation, is described in the Methods section of the manuscript.

The final CDPS model was generated using the standard workflow implemented in the Mime R package, without custom modifications to the package source code. The available object `model_rsf_ridge.rds` corresponds to the final trained model selected during the original model-development workflow and is provided here for CDPS calculation.
