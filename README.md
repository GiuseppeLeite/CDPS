# Integrated Profiling of Non-Apoptotic Cell Death Identifies SLC3A2-Mediated Ferroptosis Suppression as a Prognostic Determinant in Melanoma

This repository provides the trained Cell Death Programmed Score (CDPS) model and an example script to apply the model to the GSE91061 melanoma immunotherapy cohort.

The CDPS model was developed as described in the manuscript using genes associated with regulated non-apoptotic cell death pathways, including necroptosis, pyroptosis, and ferroptosis. The repository does not retrain the model. It provides the final trained model object and a reproducible example of CDPS application.

## Repository contents

```text
.
├── CDPS_code.R
├── res_11_list_148_rsf_ridge.rds
├── GSE91061_mime.rds
└── README.md
