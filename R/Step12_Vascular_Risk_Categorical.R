# Step12_Vascular_Risk_Categorical.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 12
source("R/lib/oef_utils.R")
oef_start_step("Step12")

# ============================================================
# Purpose:
#   Re-run the categorical VRS analysis from the original script.
#
# Analysis is intentionally kept consistent with the original script:
#   - Input: trust_vascular_risk_zscores.csv
#   - Outcome: z_score
#   - VRS = BMI_Code + HT_Code + HL_Code + DB_Code
#   - Categorical model: z_score ~ factor(VRS), levels 0-4, VRS=0 reference
# ============================================================

# ----------------------------
# 1) Paths and variable settings
# ----------------------------
InFile <- oef_input("trust_vascular_risk_zscores.csv")

RevisionRoot <- oef_output("vascular_risk_categorical")

# A NEW timestamped subfolder is created on every run.
# This guarantees that no previous output is overwritten.
RunTag <- format(Sys.time(), "%Y%m%d_%H%M%S")
OutDir <- file.path(RevisionRoot, paste0("VRS_Categorical_", RunTag))

Zvar <- "z_score"
VRS_components <- c("BMI_Code", "HT_Code", "HL_Code", "DB_Code")

ExcludeStrokeCol <- "StrokeHistory"
ExcludeHeartCol  <- "HeartDisease"

# ----------------------------
# 2) Helper functions
# ----------------------------
dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(OutDir)) {
  stop("Could not create output folder: ", OutDir)
}

format_p_text <- function(p) {
  if (is.na(p) || !is.finite(p)) return("NA")
  if (p < 0.0001) return("P < 0.0001")
  paste0("P = ", formatC(p, format = "fg", digits = 3))
}

safe_write_csv <- function(x, filename) {
  outfile <- file.path(OutDir, filename)
  if (file.exists(outfile)) {
    stop("Refusing to overwrite existing file: ", outfile)
  }
  write.csv(x, outfile, row.names = FALSE)
  message("Saved: ", outfile)
}

# ----------------------------
# 3) Load data
# ----------------------------
if (!file.exists(InFile)) {
  stop("Input file not found: ", InFile)
}

dat <- read.csv(InFile, stringsAsFactors = FALSE)

required_cols <- c(Zvar, VRS_components)
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Make sure outcome is numeric
dat[[Zvar]] <- suppressWarnings(as.numeric(dat[[Zvar]]))

cat("Original N =", nrow(dat), "\n")

# ----------------------------
# 4) Apply the SAME exclusions as the original script
# ----------------------------
dat_filt <- dat

if (ExcludeStrokeCol %in% names(dat_filt)) {
  n0 <- nrow(dat_filt)
  dat_filt <- dat_filt[
    dat_filt[[ExcludeStrokeCol]] == 0 | is.na(dat_filt[[ExcludeStrokeCol]]),
    ,
    drop = FALSE
  ]
  cat("Excluded StrokeHistory == 1:", n0 - nrow(dat_filt), "\n")
} else {
  cat("StrokeHistory column not found; no stroke exclusion applied.\n")
}

if (ExcludeHeartCol %in% names(dat_filt)) {
  n1 <- nrow(dat_filt)
  dat_filt <- dat_filt[
    dat_filt[[ExcludeHeartCol]] == 0 | is.na(dat_filt[[ExcludeHeartCol]]),
    ,
    drop = FALSE
  ]
  cat("Excluded HeartDisease == 1:", n1 - nrow(dat_filt), "\n")
} else {
  cat("HeartDisease column not found; no heart-disease exclusion applied.\n")
}

cat("N after exclusions =", nrow(dat_filt), "\n")

# ----------------------------
# 5) Recalculate VRS exactly as before
# ----------------------------
for (v in VRS_components) {
  new_name <- paste0(v, "_num")
  dat_filt[[new_name]] <- oef_binary(dat_filt[[v]], v)
  dat_filt[[new_name]][!(dat_filt[[new_name]] %in% c(0, 1))] <- NA_integer_
}

vrs_num_cols <- paste0(VRS_components, "_num")

dat_filt$VRS_complete <- complete.cases(dat_filt[, vrs_num_cols, drop = FALSE])

dat_filt$VRS_numeric <- NA_integer_
dat_filt$VRS_numeric[dat_filt$VRS_complete] <- rowSums(
  dat_filt[dat_filt$VRS_complete, vrs_num_cols, drop = FALSE]
)

# Analysis sample: complete VRS + finite z-score
dat_vrs <- dat_filt[
  dat_filt$VRS_complete &
    is.finite(dat_filt$VRS_numeric) &
    is.finite(dat_filt[[Zvar]]),
  ,
  drop = FALSE
]

dat_vrs$VRS_cat <- factor(
  as.integer(round(dat_vrs$VRS_numeric)),
  levels = 0:4
)

dat_vrs$VRS_cat <- droplevels(dat_vrs$VRS_cat)

if (nrow(dat_vrs) < 5 || nlevels(dat_vrs$VRS_cat) < 2) {
  stop("Insufficient data for categorical VRS analysis.")
}

if (!("0" %in% levels(dat_vrs$VRS_cat))) stop("VRS = 0 is required as the reference category.")

cat("Final categorical-analysis N =", nrow(dat_vrs), "\n")
cat("VRS distribution:\n")
print(table(dat_vrs$VRS_cat))

# ----------------------------
# 6) Group-wise descriptive results
#    N, mean, SD, SE, and 95% CI of mean z-score
# ----------------------------
vrs_levels <- levels(dat_vrs$VRS_cat)

group_summary_list <- lapply(vrs_levels, function(g) {
  x <- dat_vrs[[Zvar]][dat_vrs$VRS_cat == g]
  x <- x[is.finite(x)]

  n <- length(x)
  mn <- mean(x)
  s <- sd(x)
  se <- s / sqrt(n)

  if (n > 1) {
    crit <- qt(0.975, df = n - 1)
    ci_low <- mn - crit * se
    ci_high <- mn + crit * se
  } else {
    ci_low <- NA_real_
    ci_high <- NA_real_
  }

  data.frame(
    VRS = as.integer(as.character(g)),
    N = n,
    Mean_z = mn,
    SD_z = s,
    SE_z = se,
    Mean_95CI_low = ci_low,
    Mean_95CI_high = ci_high,
    stringsAsFactors = FALSE
  )
})

group_summary <- do.call(rbind, group_summary_list)
safe_write_csv(group_summary, "01_VRS_categorical_group_summary.csv")

# ----------------------------
# 7) Main categorical model
#    VRS=0 is the reference group
# ----------------------------
dat_vrs$VRS_cat <- relevel(dat_vrs$VRS_cat, ref = "0")

fit_cat <- lm(
  as.formula(paste0(Zvar, " ~ VRS_cat")),
  data = dat_vrs
)

cat("\n===== CATEGORICAL MODEL =====\n")
print(summary(fit_cat))

# Omnibus ANOVA: tests whether ANY VRS category differs
aov_cat <- anova(fit_cat)

df1 <- aov_cat["VRS_cat", "Df"]
df2 <- aov_cat["Residuals", "Df"]
F_value <- aov_cat["VRS_cat", "F value"]
P_omnibus <- aov_cat["VRS_cat", "Pr(>F)"]

omnibus_result <- data.frame(
  Model = "OEF z-score ~ categorical VRS (0-4)",
  N = nrow(dat_vrs),
  Reference = "VRS 0",
  Numerator_df = df1,
  Denominator_df = df2,
  F_value = F_value,
  P_value = P_omnibus,
  stringsAsFactors = FALSE
)

safe_write_csv(omnibus_result, "02_VRS_categorical_omnibus_ANOVA.csv")

# ----------------------------
# 8) Category-specific coefficients vs VRS=0
#    beta = mean difference in OEF z-score vs VRS 0
# ----------------------------
coef_tab <- summary(fit_cat)$coefficients
ci_tab <- confint(fit_cat, level = 0.95)

coef_rows <- grep("^VRS_cat", rownames(coef_tab))

category_coefficients <- data.frame(
  Term = rownames(coef_tab)[coef_rows],
  Contrast = paste0(
    "VRS ",
    sub("^VRS_cat", "", rownames(coef_tab)[coef_rows]),
    " vs VRS 0"
  ),
  Beta = coef_tab[coef_rows, "Estimate"],
  SE = coef_tab[coef_rows, "Std. Error"],
  t_value = coef_tab[coef_rows, "t value"],
  P_value = coef_tab[coef_rows, "Pr(>|t|)"],
  CI_95_low = ci_tab[coef_rows, 1],
  CI_95_high = ci_tab[coef_rows, 2],
  stringsAsFactors = FALSE
)

safe_write_csv(
  category_coefficients,
  "03_VRS_categorical_coefficients_vs_VRS0.csv"
)

# ----------------------------
# 9) Continuous VRS model for cross-check with primary analysis
#    This is NOT the categorical sensitivity analysis;
#    it is saved only to verify consistency with the manuscript.
# ----------------------------
fit_cont <- lm(
  as.formula(paste0(Zvar, " ~ VRS_numeric")),
  data = dat_vrs
)

cont_coef <- summary(fit_cont)$coefficients["VRS_numeric", ]
cont_ci <- confint(fit_cont, "VRS_numeric", level = 0.95)

continuous_result <- data.frame(
  Model = "OEF z-score ~ continuous VRS",
  N = nrow(dat_vrs),
  Beta_per_1_VRS = cont_coef["Estimate"],
  SE = cont_coef["Std. Error"],
  t_value = cont_coef["t value"],
  P_value = cont_coef["Pr(>|t|)"],
  CI_95_low = cont_ci[1],
  CI_95_high = cont_ci[2],
  stringsAsFactors = FALSE
)

safe_write_csv(
  continuous_result,
  "04_VRS_continuous_model_crosscheck.csv"
)

# ----------------------------
# 10) Model fit summary
# ----------------------------
s_cat <- summary(fit_cat)

model_fit <- data.frame(
  Model = "Categorical VRS 0-4",
  N = nrow(dat_vrs),
  R_squared = s_cat$r.squared,
  Adjusted_R_squared = s_cat$adj.r.squared,
  Residual_SE = s_cat$sigma,
  stringsAsFactors = FALSE
)

safe_write_csv(model_fit, "05_VRS_categorical_model_fit.csv")

# ----------------------------
# 11) Human-readable text report
#     Convenient for analysis report / manuscript Results
# ----------------------------
mean_text <- paste(
  sprintf(
    "VRS %d: %.3f +/- %.3f (n=%d)",
    group_summary$VRS,
    group_summary$Mean_z,
    group_summary$SD_z,
    group_summary$N
  ),
  collapse = "; "
)

contrast_text <- paste(
  sprintf(
    "%s: beta=%.3f, 95%% CI [%.3f, %.3f], %s",
    category_coefficients$Contrast,
    category_coefficients$Beta,
    category_coefficients$CI_95_low,
    category_coefficients$CI_95_high,
    vapply(category_coefficients$P_value, format_p_text, character(1))
  ),
  collapse = "; "
)

report_lines <- c(
  "VRS CATEGORICAL SENSITIVITY ANALYSIS",
  "====================================",
  "",
  paste0("Input: ", InFile),
  paste0("Output folder: ", OutDir),
  paste0("Final N: ", nrow(dat_vrs)),
  "",
  "Omnibus categorical test:",
  sprintf(
    "F(%d, %d) = %.3f, %s",
    as.integer(df1),
    as.integer(df2),
    F_value,
    format_p_text(P_omnibus)
  ),
  "",
  "Group-wise OEF z-score mean +/- SD:",
  mean_text,
  "",
  "Category-specific differences relative to VRS=0:",
  contrast_text,
  "",
  "Continuous VRS cross-check:",
  sprintf(
    "beta = %.3f, SE = %.3f, 95%% CI [%.3f, %.3f], %s",
    continuous_result$Beta_per_1_VRS,
    continuous_result$SE,
    continuous_result$CI_95_low,
    continuous_result$CI_95_high,
    format_p_text(continuous_result$P_value)
  ),
  "",
  "Suggested Results wording:",
  sprintf(
    paste0(
      "When VRS was modeled categorically rather than continuously, ",
      "OEF z-scores remained significantly different across vascular-risk ",
      "burden categories (F(%d, %d) = %.2f, %s)."
    ),
    as.integer(df1),
    as.integer(df2),
    F_value,
    format_p_text(P_omnibus)
  )
)

report_file <- file.path(OutDir, "06_VRS_categorical_analysis_report.txt")
if (file.exists(report_file)) {
  stop("Refusing to overwrite existing file: ", report_file)
}
writeLines(report_lines, con = report_file)

cat("\n============================================================\n")
cat("Analysis completed successfully.\n")
cat("NO previous files were modified or overwritten.\n")
cat("All new outputs are in:\n", OutDir, "\n", sep = "")
cat("============================================================\n")
