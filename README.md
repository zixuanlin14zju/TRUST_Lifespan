# TRUST Lifespan: OEF analysis code

Analysis scripts for lifespan normative modeling of oxygen extraction fraction (OEF), clinical applications, reproducibility, and physiological sensitivity analyses. The release retains the original analysis sequence as **Step0–Step7**, adds the supplied supplementary analyses as **Step8–Step20**, and numbers the previous resampling scripts **Step21–Step26**. 

## Scope and release status

**This repository distributes code, not participant data.** No original or sampled study records, clinical workbooks, fitted models, historical plots, or R session histories are included. Data access remains subject to the contributing institutions' approvals and governance policies.

## Analysis map

| Step | Analysis | Inputs / dependency |
|---|---|---|
| 0 | Family and cubic spline-complexity comparison | Normative HC table |
| 1 | Main BCTo lifespan model, centiles, variability, CV and z-scores | Normative HC table |
| 2 | Vascular risk, individual factors, multivariable factors, categorical sensitivity | Approved HC z-score/vascular-risk table |
| 3 | Disease deviation patterns | HC and patient tables; refits the supplied fixed-complexity BCTo model |
| 4 | Tumor phenotype analyses | Processed tumor table with z-scores |
| 5 | Exploratory cognition/APOE analyses | Processed neurodegenerative table with z-scores and APOE_Code |
| 6 | Caffeine follow-up reproducibility | Paired-day OEF table |
| 7 | Traveling OEF reproducibility | Original paired-site OEF table |
| 8 | Four Hct normative sensitivity models | Precomputed HC and patient blood-sensitivity workbooks |
| 9 | Hybrid pediatric/adult Ya sensitivity | Precomputed HC blood-sensitivity workbook |
| 10 | Exploratory Hct sensitivity of clinical deviations | HC and patient blood-sensitivity workbooks |
| 11 | Additional within-disease age regression | HC and patient tables |
| 12 | Categorical VRS contrasts and report | HC z-score/vascular-risk table |
| 13 | Community aging: VRS, CBF, OEF and CMRO2 | Community-aging paired-CBF workbook |
| 14 | Pediatric OSA: CBF/CMRO2 and CBF-adjusted OEF | Pediatric OSA/control workbook |
| 15 | Cognitive aging: diagnosis, cognition, CBF, VRS, CMRO2 | Cognitive-aging workbook |
| 16 | Adults at Site 1 versus Sites 8/12 | HC workbook; common ages 23–33 years |
| 17 | Early-life QC and mu-spline complexity | HC workbook with dR2 |
| 18 | Direct venous R2 distribution by site | HC workbook with T2 in milliseconds |
| 19 | Hct/Ya trajectory and derivative presentation | Outputs from Steps 8 and 9 |
| 20 | Reference Hct/Ya assumption curves | Constants in the supplied plotting script; no participant data |
| 21 | Age-balanced resampling | Normative HC table |
| 22 | Bootstrap | Normative HC table |
| 23 | Split-half stability | Normative HC table |
| 24 | Leave-one-site-out stability | Normative HC table |
| 25 | Stricter-QC sensitivity | Independently quality-controlled HC table |
| 26 | Numerical concordance of sensitivity curves | Normative HC table and outputs from Steps 21–25 |
| 27 | Paired OEF and optional R2 site comparisons | Same input as Step7 |

TRUST quantification is outside this R analysis repository and will be supplied separately in MATLAB. The R scripts start from quantified physiological variables. 
