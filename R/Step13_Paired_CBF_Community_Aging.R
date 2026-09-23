# Step13_Paired_CBF_Community_Aging.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 13
source("R/lib/oef_utils.R")
oef_start_step("Step13")

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

suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
})

# ============================================================
# 2. PATHS
# ============================================================

DataDir <-
  "data"

# Use the explicitly named cohort input.
InputFile <- oef_input("trust_community_aging_cbf.xlsx")
if (!file.exists(InputFile)) stop("Required private input is missing: ", InputFile)
OutputDir <- file.path(oef_output(), "paired_cbf_community")
dir.create(OutputDir, recursive = TRUE, showWarnings = FALSE)

cat("\n============================================================\n")
cat("INPUT FILE\n")
cat("============================================================\n")

cat(InputFile, "\n")

# ============================================================
# 3. READ DATA
# ============================================================

dat <- read_excel(
  InputFile,
  sheet = "Sheet1"
)

cat("\nOriginal N =", nrow(dat), "\n")

# ============================================================
# 4. VARIABLE CHECK
# ============================================================

RequiredVariables <- c(
  "Age",
  "Sex",
  "OEF zscore",
  "CBF",
  "CMRO2",
  "BMI_Code",
  "HT_Code",
  "HL_Code",
  "DB_Code"
)

MissingVariables <- setdiff(
  RequiredVariables,
  names(dat)
)

if (length(MissingVariables) > 0) {

  stop(
    "Missing required variable(s): ",
    paste(
      MissingVariables,
      collapse = ", "
    )
  )
}

# ============================================================
# 5. RENAME VARIABLES
# ============================================================

dat$OEFz <- suppressWarnings(
  as.numeric(dat[["OEF zscore"]])
)

dat$CBF <- suppressWarnings(
  as.numeric(dat$CBF)
)

dat$CMRO2 <- suppressWarnings(
  as.numeric(dat$CMRO2)
)

dat$Age <- suppressWarnings(
  as.numeric(dat$Age)
)

# ============================================================
# 6. SEX
#
# Current cohort coding:
# 1 = Male
# 2 = Female
# ============================================================

dat$SexF <- factor(
  dat$Sex,
  levels = c(1, 2),
  labels = c(
    "Male",
    "Female"
  )
)

cat("\nSex distribution:\n")
print(
  table(
    dat$SexF,
    useNA = "ifany"
  )
)

# ============================================================
# 7. CONSTRUCT FOUR-COMPONENT VRS
#
# Current TRUST lifespan manuscript definition:
#
# BMI_Code
# HT_Code
# HL_Code
# DB_Code
#
# Each:
# 0 = absent
# 1 = present
#
# VRS range = 0-4
#
# IMPORTANT:
# If any component is missing, VRS is missing.
# ============================================================

VRS_components <- c(
  "BMI_Code",
  "HT_Code",
  "HL_Code",
  "DB_Code"
)

# Convert explicitly to numeric
for (v in VRS_components) {

  dat[[v]] <- suppressWarnings(
    as.numeric(dat[[v]])
  )
}

dat$VRS <- rowSums(
  dat[, VRS_components],
  na.rm = FALSE
)

# Check whether VRS contains unexpected values
if (
  any(
    !is.na(dat$VRS) &
    !(dat$VRS %in% 0:4)
  )
) {

  warning(
    "Some VRS values are outside the expected 0-4 range."
  )
}

cat("\n============================================================\n")
cat("VRS DISTRIBUTION\n")
cat("============================================================\n")

print(
  table(
    dat$VRS,
    useNA = "ifany"
  )
)

tidy_lm <- function(
    model,
    model_name
) {

  sm <- summary(model)$coefficients

  ci <- confint(
    model,
    level = 0.95
  )

  out <- data.frame(

    Model = model_name,

    Term =
      rownames(sm),

    Estimate =
      sm[, "Estimate"],

    SE =
      sm[, "Std. Error"],

    CI_low =
      ci[, 1],

    CI_high =
      ci[, 2],

    t_value =
      sm[, "t value"],

    df =
      df.residual(model),

    P_value =
      sm[, "Pr(>|t|)"],

    N =
      nobs(model),

    R2 =
      summary(model)$r.squared,

    Adjusted_R2 =
      summary(model)$adj.r.squared,

    stringsAsFactors = FALSE
  )

  rownames(out) <- NULL

  out
}

extract_term <- function(
    model,
    term_name,
    description
) {

  sm <- summary(model)$coefficients
  ci <- confint(model)

  if (!(term_name %in% rownames(sm))) {
    stop(
      "Cannot find term: ",
      term_name
    )
  }

  data.frame(

    Analysis =
      description,

    N =
      nobs(model),

    Estimate =
      sm[term_name, "Estimate"],

    SE =
      sm[term_name, "Std. Error"],

    CI_low =
      ci[term_name, 1],

    CI_high =
      ci[term_name, 2],

    t_value =
      sm[term_name, "t value"],

    df =
      df.residual(model),

    P_value =
      sm[term_name, "Pr(>|t|)"],

    stringsAsFactors = FALSE
  )
}

fit_CBF_VRS <- lm(

  CBF ~
    VRS +
    Age +
    SexF,

  data = dat,

  na.action = na.omit
)

fit_CMRO2_VRS <- lm(

  CMRO2 ~
    VRS +
    Age +
    SexF,

  data = dat,

  na.action = na.omit
)

fit_OEFz_VRS_CBF <- lm(

  OEFz ~
    VRS +
    CBF +
    Age +
    SexF,

  data = dat,

  na.action = na.omit
)

dat_same <- dat[
  complete.cases(
    dat[, c(
      "VRS",
      "OEFz",
      "CBF",
      "Age",
      "SexF"
    )]
  ),
]

fit_OEFz_base <- lm(

  OEFz ~
    CBF +
    Age +
    SexF,

  data = dat_same
)

# Reuse the retained full model rather than fit the same design twice.
# dat_same is exactly its complete-case sample; missing CMRO2 is NOT a filter.
fit_OEFz_addVRS <- fit_OEFz_VRS_CBF
stopifnot(nobs(fit_OEFz_base) == nobs(fit_OEFz_addVRS))
stopifnot(isTRUE(all.equal(
  model.frame(fit_OEFz_base),
  model.frame(fit_OEFz_addVRS)[, names(model.frame(fit_OEFz_base)), drop = FALSE],
  check.attributes = FALSE
)))

incremental_anova <-
  anova(
    fit_OEFz_base,
    fit_OEFz_addVRS
  )

incremental_result <- data.frame(

  Base_Model =
    "OEFz ~ CBF + Age + Sex",

  Added_Model =
    "OEFz ~ CBF + Age + Sex + VRS",

  N =
    nobs(fit_OEFz_addVRS),

  Base_R2 =
    summary(fit_OEFz_base)$r.squared,

  Added_R2 =
    summary(fit_OEFz_addVRS)$r.squared,

  Delta_R2 =
    summary(fit_OEFz_addVRS)$r.squared -
    summary(fit_OEFz_base)$r.squared,

  Base_Adjusted_R2 =
    summary(fit_OEFz_base)$adj.r.squared,

  Added_Adjusted_R2 =
    summary(fit_OEFz_addVRS)$adj.r.squared,

  F_change =
    incremental_anova$F[2],

  P_change =
    incremental_anova$`Pr(>F)`[2],

  stringsAsFactors = FALSE
)

# --- Only the three reported coefficients; Delta R2 is in Incremental_Test ---
main_results <- rbind(
  extract_term(fit_CBF_VRS, "VRS", "VRS -> CBF, adjusted for age and sex"),
  extract_term(fit_OEFz_VRS_CBF, "VRS",
               "VRS -> OEF z-score, adjusted for CBF, age, and sex"),
  extract_term(fit_CMRO2_VRS, "VRS", "VRS -> CMRO2, adjusted for age and sex")
)
all_models <- rbind(
  tidy_lm(fit_CBF_VRS, "CBF ~ VRS + Age + Sex"),
  tidy_lm(fit_OEFz_VRS_CBF, "OEFz ~ VRS + CBF + Age + Sex"),
  tidy_lm(fit_CMRO2_VRS, "CMRO2 ~ VRS + Age + Sex"),
  tidy_lm(fit_OEFz_base, "OEFz ~ CBF + Age + Sex [same complete-case sample]")
)
sample_summary <- data.frame(
  Analysis = c("Input rows", "Complete VRS", "CBF model", "OEFz + CBF model",
               "CMRO2 model", "Delta R2: common complete-case sample"),
  N = c(nrow(dat), sum(!is.na(dat$VRS)), nobs(fit_CBF_VRS),
        nobs(fit_OEFz_VRS_CBF), nobs(fit_CMRO2_VRS), nrow(dat_same))
)
cat("\n MAIN RESULTS\n")
print(main_results, row.names = FALSE)
cat("\nINCREMENTAL VRS CONTRIBUTION\n")
print(incremental_result, row.names = FALSE)

# Keep analysis outputs separate from historical exploratory workbooks.
TimeStamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OutputFile <- file.path(OutputDir,
  paste0("Community_VRS_CBF_OEF_", TimeStamp, ".xlsx"))
wb <- createWorkbook()
output_sheets <- list(
  Main_Results = main_results,
  Incremental_Test = incremental_result,
  All_Models = all_models,
  Sample_Summary = sample_summary
)
header_style <- createStyle(textDecoration = "bold", border = "Bottom")
for (sheet_name in names(output_sheets)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, output_sheets[[sheet_name]])
  addStyle(wb, sheet_name, header_style, rows = 1,
           cols = seq_len(ncol(output_sheets[[sheet_name]])), gridExpand = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(output_sheets[[sheet_name]])),
               widths = 18)
  setColWidths(wb, sheet_name, cols = 1, widths = 58)
  freezePane(wb, sheet_name, firstRow = TRUE)
}
saveWorkbook(wb, OutputFile, overwrite = FALSE)
cat("\nResults saved to:\n", OutputFile, "\n")
