# Step14_Paired_CBF_Pediatric_OSA.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 14
source("R/lib/oef_utils.R")
oef_start_step("Step14")


# ============================================================
# 0. CLEAR WORKSPACE
# ============================================================

# Run each Step in a separate R session; no workspace deletion.
gc()

options(stringsAsFactors = FALSE)

# ============================================================
# 1. PACKAGES
# ============================================================

required_packages <- c(
  "readxl",
  "openxlsx"
)

missing_packages <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "), ". Run Rscript setup.R --install explicitly.")
}

library(readxl)
library(openxlsx)

# ============================================================
# 2. PATH
# ============================================================

DataDir <- oef_input()

# Use the explicitly named cohort input.
InputFile <- oef_input("trust_pediatric_osa_cbf.xlsx")
if (!file.exists(InputFile)) stop("Required private input is missing: ", InputFile)
OutputDir <- file.path(oef_output(), "paired_cbf_osa")
dir.create(OutputDir, recursive = TRUE, showWarnings = FALSE)

cat("\nUsing input file:\n")
cat(InputFile, "\n\n")

# ============================================================
# 3. READ DATA
# ============================================================

dat <- read_excel(
  InputFile,
  sheet = 1
)

cat("Original N =", nrow(dat), "\n\n")

cat("Variables in input file:\n")
print(names(dat))

# ============================================================
# 4. STANDARDIZE VARIABLE NAMES
# ============================================================

# Sex: original coding M = 0, F = 1
if ("Sex(M0F1)" %in% names(dat)) {
  dat$Sex <- dat[["Sex(M0F1)"]]
} else if (!("Sex" %in% names(dat))) {
  stop("Cannot find sex variable: expected Sex(M0F1) or Sex.")
}

# OEF z-score
if ("OEFzscore" %in% names(dat)) {
  dat$OEFz <- dat$OEFzscore
} else if ("OEF zscore" %in% names(dat)) {
  dat$OEFz <- dat[["OEF zscore"]]
} else if (!("OEFz" %in% names(dat))) {
  stop("Cannot find OEF z-score variable.")
}

# Diagnosis
if (!("Diagnosis" %in% names(dat))) {
  stop("Cannot find Diagnosis variable.")
}

# CMRO2
if (!("CMRO2" %in% names(dat))) {
  stop("Cannot find CMRO2 variable.")
}

# ============================================================
# 5. CREATE OSA INDICATOR
#
# Normal = 0
# OSA    = 1
# ============================================================

dat$Diagnosis_clean <- toupper(
  trimws(as.character(dat$Diagnosis))
)

dat$OSA <- ifelse(
  dat$Diagnosis_clean == "OSA",
  1,
  ifelse(
    dat$Diagnosis_clean == "NORMAL",
    0,
    NA
  )
)

cat("\nDiagnosis counts:\n")
print(table(
  dat$Diagnosis_clean,
  useNA = "ifany"
))

# ============================================================
# 6. ANALYSIS DATASETS
#
# Use model-specific complete cases so missing CMRO2 does not
# unnecessarily remove subjects from the CBF/OEF analyses.
# ============================================================

ana_CBF <- dat[
  complete.cases(
    dat[, c(
      "CBF",
      "Age",
      "Sex",
      "OSA"
    )]
  ),
]

ana_CMRO2 <- dat[
  complete.cases(
    dat[, c(
      "CMRO2",
      "Age",
      "Sex",
      "OSA"
    )]
  ),
]

ana_OEF <- dat[
  complete.cases(
    dat[, c(
      "OEFz",
      "CBF",
      "Age",
      "Sex",
      "OSA"
    )]
  ),
]

cat("\n============================================================\n")
cat("SAMPLE SIZE\n")
cat("============================================================\n")

cat(
  "CBF model: N =",
  nrow(ana_CBF),
  "; Normal =",
  sum(ana_CBF$OSA == 0),
  "; OSA =",
  sum(ana_CBF$OSA == 1),
  "\n"
)

cat(
  "CMRO2 model: N =",
  nrow(ana_CMRO2),
  "; Normal =",
  sum(ana_CMRO2$OSA == 0),
  "; OSA =",
  sum(ana_CMRO2$OSA == 1),
  "\n"
)

cat(
  "OEF model: N =",
  nrow(ana_OEF),
  "; Normal =",
  sum(ana_OEF$OSA == 0),
  "; OSA =",
  sum(ana_OEF$OSA == 1),
  "\n"
)

# ============================================================
# 7. HELPER FUNCTION FOR MODEL RESULT
# ============================================================

extract_lm_term <- function(
    model,
    term,
    description
) {

  sm <- summary(model)$coefficients
  ci <- confint(model)

  if (!(term %in% rownames(sm))) {
    stop("Cannot find term ", term)
  }

  data.frame(
    Analysis = description,
    N = nobs(model),
    Estimate = sm[term, "Estimate"],
    SE = sm[term, "Std. Error"],
    CI_low = ci[term, 1],
    CI_high = ci[term, 2],
    t_value = sm[term, "t value"],
    df = df.residual(model),
    P_value = sm[term, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}

tidy_model <- function(
    model,
    model_name
) {

  sm <- summary(model)$coefficients
  ci <- confint(model)

  out <- data.frame(
    Model = model_name,
    Term = rownames(sm),
    Estimate = sm[, "Estimate"],
    SE = sm[, "Std. Error"],
    CI_low = ci[, 1],
    CI_high = ci[, 2],
    t_value = sm[, "t value"],
    P_value = sm[, "Pr(>|t|)"],
    N = nobs(model),
    stringsAsFactors = FALSE
  )

  rownames(out) <- NULL

  out
}

# ============================================================
# 8. ANALYSIS 1
# CBF ~ OSA + AGE + SEX
#
# Question:
# Is resting cerebral perfusion altered in OSA?
# ============================================================

fit_CBF <- lm(
  CBF ~ OSA + Age + Sex,
  data = ana_CBF
)

cat("\n============================================================\n")
cat("ANALYSIS 1: CBF ~ OSA + Age + Sex\n")
cat("============================================================\n\n")

print(summary(fit_CBF))

result_CBF <- extract_lm_term(
  fit_CBF,
  "OSA",
  "CBF: OSA vs Normal, adjusted for age and sex"
)

cat("\nAdjusted OSA effect on CBF:\n")

cat(
  sprintf(
    paste0(
      "beta = %.3f, ",
      "95%% CI [%.3f, %.3f], ",
      "P = %.5f\n"
    ),
    result_CBF$Estimate,
    result_CBF$CI_low,
    result_CBF$CI_high,
    result_CBF$P_value
  )
)

# ============================================================
# 9. ANALYSIS 2
# CMRO2 ~ OSA + AGE + SEX
#
# Question:
# Is cerebral oxygen metabolic rate altered in OSA?
# ============================================================

fit_CMRO2 <- lm(
  CMRO2 ~ OSA + Age + Sex,
  data = ana_CMRO2
)

cat("\n============================================================\n")
cat("ANALYSIS 2: CMRO2 ~ OSA + Age + Sex\n")
cat("============================================================\n\n")

print(summary(fit_CMRO2))

result_CMRO2 <- extract_lm_term(
  fit_CMRO2,
  "OSA",
  "CMRO2: OSA vs Normal, adjusted for age and sex"
)

cat("\nAdjusted OSA effect on CMRO2:\n")

cat(
  sprintf(
    paste0(
      "beta = %.3f, ",
      "95%% CI [%.3f, %.3f], ",
      "P = %.5f\n"
    ),
    result_CMRO2$Estimate,
    result_CMRO2$CI_low,
    result_CMRO2$CI_high,
    result_CMRO2$P_value
  )
)

# ============================================================
# 10. NORMAL-REFERENCE CENTERING FOR OEF MODEL
#
# Center CBF, age, and sex at the mean values of NORMAL controls.
#
# This makes:
#
# Intercept:
# adjusted mean OEF z in normal controls
#
# Intercept + OSA coefficient:
# adjusted mean OEF z in OSA
#
# at the reference CBF / age / sex distribution.
# ============================================================

normal_OEF <- ana_OEF[
  ana_OEF$OSA == 0,
]

ref_CBF <- mean(
  normal_OEF$CBF,
  na.rm = TRUE
)

ref_Age <- mean(
  normal_OEF$Age,
  na.rm = TRUE
)

ref_Sex <- mean(
  normal_OEF$Sex,
  na.rm = TRUE
)

ana_OEF$CBF_c <-
  ana_OEF$CBF - ref_CBF

ana_OEF$Age_c <-
  ana_OEF$Age - ref_Age

ana_OEF$Sex_c <-
  ana_OEF$Sex - ref_Sex

cat("\n============================================================\n")
cat("NORMAL REFERENCE VALUES FOR OEF MODEL\n")
cat("============================================================\n")

cat("Mean CBF =", ref_CBF, "\n")
cat("Mean Age =", ref_Age, "\n")
cat(
  "Mean Sex =",
  ref_Sex,
  "(female proportion because M=0, F=1)\n"
)

# ============================================================
# 11. ANALYSIS 3
# OEF z ~ OSA + CBF + AGE + SEX
#
# Main supplementary test:
# Does OEF provide information beyond CBF?
# ============================================================

fit_OEF <- lm(
  OEFz ~
    OSA +
    CBF_c +
    Age_c +
    Sex_c,
  data = ana_OEF
)

cat("\n============================================================\n")
cat("ANALYSIS 3: OEF z ~ OSA + CBF + Age + Sex\n")
cat("============================================================\n\n")

print(summary(fit_OEF))

result_OEF_OSA <- extract_lm_term(
  fit_OEF,
  "OSA",
  paste0(
    "OEF z: OSA vs Normal, ",
    "adjusted for CBF, age, and sex"
  )
)

cat("\nAdjusted OSA effect on OEF z-score:\n")

cat(
  sprintf(
    paste0(
      "beta = %.3f, ",
      "95%% CI [%.3f, %.3f], ",
      "P = %.5f\n"
    ),
    result_OEF_OSA$Estimate,
    result_OEF_OSA$CI_low,
    result_OEF_OSA$CI_high,
    result_OEF_OSA$P_value
  )
)

# ============================================================
# 12. RAW GROUP SUMMARY
# ============================================================

group_summary_function <- function(
    data,
    variable
) {

  # Check variable exists
  if (!(variable %in% names(data))) {
    stop("Variable not found: ", variable)
  }

  # Extract as a plain numeric vector
  x_all <- suppressWarnings(
    as.numeric(data[[variable]])
  )

  osa_all <- data$OSA

  # Keep rows with both outcome and group available
  keep <- is.finite(x_all) & !is.na(osa_all)

  x_all <- x_all[keep]
  osa_all <- osa_all[keep]

  out <- do.call(
    rbind,
    lapply(
      c(0, 1),
      function(g) {

        # IMPORTANT:
        # extract numeric vector, not one-column tibble
        x <- x_all[osa_all == g]

        data.frame(
          Variable = variable,
          Group = ifelse(
            g == 0,
            "Normal",
            "OSA"
          ),
          N = length(x),
          Mean = mean(x, na.rm = TRUE),
          SD = sd(x, na.rm = TRUE),
          Median = median(x, na.rm = TRUE),
          Min = min(x, na.rm = TRUE),
          Max = max(x, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(out) <- NULL

  return(out)
}

# ------------------------------------------------------------
# Generate group summaries
# ------------------------------------------------------------

group_summary <- rbind(

  group_summary_function(
    dat,
    "CBF"
  ),

  group_summary_function(
    dat,
    "CMRO2"
  ),

  group_summary_function(
    dat,
    "OEFz"
  )
)

cat("\n============================================================\n")
cat("GROUP SUMMARY\n")
cat("============================================================\n\n")

print(group_summary)

# ============================================================
# 13. FULL MODEL COEFFICIENT TABLES
# ============================================================

model_CBF_table <- tidy_model(
  fit_CBF,
  "CBF ~ OSA + Age + Sex"
)

model_CMRO2_table <- tidy_model(
  fit_CMRO2,
  "CMRO2 ~ OSA + Age + Sex"
)

model_OEF_table <- tidy_model(
  fit_OEF,
  "OEFz ~ OSA + CBF + Age + Sex"
)

# ============================================================
# 14. FINAL MAIN RESULTS
# ============================================================

main_results <- rbind(

  result_CBF,

  result_CMRO2,

  result_OEF_OSA
)

cat("\n============================================================\n")
cat("FINAL MAIN RESULTS\n")
cat("============================================================\n\n")

print(main_results)

# ============================================================
# 15. REFERENCE VALUES
# ============================================================

reference_values <- data.frame(

  Variable = c(
    "CBF",
    "Age",
    "Sex (female proportion)"
  ),

  Normal_reference_mean = c(
    ref_CBF,
    ref_Age,
    ref_Sex
  ),

  stringsAsFactors = FALSE
)

# ============================================================
# 16. SAMPLE SUMMARY
# ============================================================

sample_summary <- data.frame(

  Analysis = c(
    "CBF model",
    "CMRO2 model",
    "OEF z / CBF model"
  ),

  Total_N = c(
    nrow(ana_CBF),
    nrow(ana_CMRO2),
    nrow(ana_OEF)
  ),

  Normal_N = c(
    sum(ana_CBF$OSA == 0),
    sum(ana_CMRO2$OSA == 0),
    sum(ana_OEF$OSA == 0)
  ),

  OSA_N = c(
    sum(ana_CBF$OSA == 1),
    sum(ana_CMRO2$OSA == 1),
    sum(ana_OEF$OSA == 1)
  ),

  stringsAsFactors = FALSE
)

# ============================================================
# 17. SAVE RESULTS TO EXCEL
# ============================================================

TimeStamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

OutputFile <- file.path(
  OutputDir,
  paste0(
    "OSA_CBF_CMRO2_OEF_",
    TimeStamp,
    ".xlsx"
  )
)

wb <- createWorkbook()

# ------------------------------------------------------------
# Main results
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Main_Results"
)

writeData(
  wb,
  "Main_Results",
  main_results
)

# ------------------------------------------------------------
# Group descriptive summary
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Group_Summary"
)

writeData(
  wb,
  "Group_Summary",
  group_summary
)

# ------------------------------------------------------------
# CBF model
# ------------------------------------------------------------

addWorksheet(
  wb,
  "CBF_Model"
)

writeData(
  wb,
  "CBF_Model",
  model_CBF_table
)

# ------------------------------------------------------------
# CMRO2 model
# ------------------------------------------------------------

addWorksheet(
  wb,
  "CMRO2_Model"
)

writeData(
  wb,
  "CMRO2_Model",
  model_CMRO2_table
)

# ------------------------------------------------------------
# OEF model
# ------------------------------------------------------------

addWorksheet(
  wb,
  "OEFz_Model"
)

writeData(
  wb,
  "OEFz_Model",
  model_OEF_table
)

# ------------------------------------------------------------
# Reference values
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Reference_Values"
)

writeData(
  wb,
  "Reference_Values",
  reference_values
)

# ------------------------------------------------------------
# Sample summary
# ------------------------------------------------------------

addWorksheet(
  wb,
  "Sample_Summary"
)

writeData(
  wb,
  "Sample_Summary",
  sample_summary
)

# ============================================================
# 18. EXCEL FORMATTING
# ============================================================

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom"
)

for (sheet_name in names(wb)) {

  sheet_data <- readWorkbook(
    wb,
    sheet = sheet_name
  )

  if (ncol(sheet_data) > 0) {

    addStyle(
      wb,
      sheet = sheet_name,
      style = header_style,
      rows = 1,
      cols = 1:ncol(sheet_data),
      gridExpand = TRUE
    )

    setColWidths(
      wb,
      sheet = sheet_name,
      cols = 1:ncol(sheet_data),
      widths = "auto"
    )

    freezePane(
      wb,
      sheet = sheet_name,
      firstRow = TRUE
    )
  }
}

# ============================================================
# 19. SAVE
# ============================================================

saveWorkbook(
  wb,
  OutputFile,
  overwrite = FALSE
)

cat("\n============================================================\n")
cat("RESULTS SAVED SUCCESSFULLY\n")
cat("============================================================\n")

cat(OutputFile, "\n")
