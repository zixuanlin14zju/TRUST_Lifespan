# TRUST Lifespan OEF Analysis Code

This repository contains R scripts used for lifespan normative modeling of global OEF measured with TRUST MRI, vascular-risk analyses, disease-deviation analyses, clinical phenotype analyses, reproducibility plots, and sensitivity analyses.

## Folder structure

```text
R/        Analysis and plotting scripts
data/     Input data tables
outputs/  Generated results and figures
```

## Data files (not uploaded yet, please contact the corresponding author)
| File | Used by |
|---|---|
| `data/trust_lifespan_hc.csv` | Main normative modeling and sensitivity analyses |
| `data/trust_lifespan_hc_stricterQC.csv` | Stricter-QC sensitivity analysis |
| `data/trust_lifespan_hc_zscores.csv` | Disease-deviation analysis |
| `data/trust_disease_patients.csv` | Disease-deviation analysis |
| `data/trust_vascular_risk_zscores.csv` | Vascular-risk analysis |
| `data/trust_tumor_phenotypes.xlsx` | Tumor clinical phenotype analysis |
| `data/trust_neurodegenerative_phenotypes.xlsx` | APOE and cognitive phenotype analysis |
| `data/trust_caffeine_followup.xlsx` | Caffeine Followup analysis |
| `data/trust_traveling_subjects.xlsx` | Travelling study analysis |

## How to run

Open R/RStudio with the repository root as the working directory. For example:

```r
setwd("path/to/TRUST_Lifespan_GitHub_Release")
```

Then run the scripts in the `R/` folder as needed. A typical order is:

```r
source("R/Step0_OEF_GAMLSS_CompareFamily.R")
source("R/Step1_OEF_GAMLSS_zscore_NoRefSite_BCTo_Sex.R")
source("R/Step2_OEF_GAMLSS_VRS.R")
source("R/Step3_OEF_GAMLSS_Disease_Pattern_NoRefSite.R")
...
```

Sensitivity analyses can be run independently:

```r
source("R/Sensitivity1_BalancedResampling_NoRefSite.R")
source("R/Sensitivity2_Bootstrap.R")
source("R/Sensitivity3_SplitinHalf_NoRefSite.R")
source("R/Sensitivity4_LOSO.R")
source("R/Sensitivity5_StricterQC.R")
source("R/Sensitivity6_ComparewithMain.R")
```

## Notes

- Paths are relative to the repository root.
- Output files are written under `outputs/`.
