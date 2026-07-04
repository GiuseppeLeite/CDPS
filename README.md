# Integrated Profiling of Non-Apoptotic Cell Death Identifies SLC3A2-Mediated Ferroptosis Suppression as a Prognostic Determinant in Melanoma

This repository provides the trained Cell Death Programmed Score (CDPS) model and an example script to apply the model to the GSE91061 melanoma immunotherapy cohort.

The CDPS model was developed as described in the manuscript using genes associated with regulated non-apoptotic cell death pathways, including necroptosis, pyroptosis, and ferroptosis. This repository does not retrain the model. It provides the final trained model object and a reproducible example of CDPS application.

## Repository contents

```text
.
├── CDPS_code.R
├── res_11_list_148_rsf_ridge.rds
├── GSE91061_mime.rds
└── README.md
```

## Files

### `CDPS_code.R`

R script used to apply the trained CDPS model to the example GSE91061 dataset and generate the paired pre-treatment versus on-treatment plot stratified by response group.

### `res_11_list_148_rsf_ridge.rds`

Trained CDPS model object.

The model object contains the fitted ridge model and the selected lambda value used for prediction:

```r
model_rsf_ridge$glmnet.fit
model_rsf_ridge$lambda.min
```

The genes used by the model are extracted from:

```r
rownames(model_rsf_ridge$glmnet.fit$beta)
```

### `GSE91061_mime.rds`

Example input dataset used to apply the CDPS model.

The dataset must be readable with `readRDS()` and convertible to a `data.frame`.

Rows must represent samples. Columns must include the 7 genes used by the CDPS model and the following annotation columns:

```text
Patient
Group
Response
```

Expected coding:

```text
Group:
  Pre = pre-treatment sample
  On  = on-treatment sample

Response:
  PRCR = responder
  other values = non-responder
```

## Input format

The input dataset should have samples in rows and genes/metadata in columns.

Example structure:

```text
Patient   Group   Response   Gene1   Gene2   Gene3   Gene4   Gene5   Gene6   Gene7
P001      Pre     PRCR       ...     ...     ...     ...     ...     ...     ...
P001      On      PRCR       ...     ...     ...     ...     ...     ...     ...
P002      Pre     PD         ...     ...     ...     ...     ...     ...     ...
P002      On      PD         ...     ...     ...     ...     ...     ...     ...
```

All 7 CDPS genes must be present as columns. The script stops if any model gene is missing, if expression values are missing or non-numeric, or if any model gene has zero variance after scaling.

## CDPS calculation

Expression values from the 7 model genes are extracted in the exact order used by the trained model. Gene expression values are Z-score standardized within the dataset before prediction.

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

The analysis requires the following R packages:

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

Input dataset with an additional column:

```text
CDPS
```

### `GSE91061_PRE_ON_ByResponse_CDPS.pdf`

Paired pre-treatment versus on-treatment CDPS plot stratified by response group.

## Notes

This repository is intended to document and reproduce the application of the trained CDPS model to GSE91061.

The model-development strategy, training cohorts, validation cohorts, feature-selection procedure, and performance evaluation are described in the Methods section of the manuscript.
