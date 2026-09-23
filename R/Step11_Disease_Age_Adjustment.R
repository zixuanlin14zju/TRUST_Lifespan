# Step11_Disease_Age_Adjustment.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 11
source("R/lib/oef_utils.R")
oef_start_step("Step11")

# ============================================================
# OEF normative model -> patient Z-scores -> age-adjusted
# disease regression sensitivity analysis
#
# Purpose:
#   Sensitivity analysis for Fig. 2f.
#   The primary OEF Z-score is already age/sex standardized by
#   the GAMLSS normative model. Here, chronological age is added
#   again as an explicit covariate within each disease group.
#
# Key point:
#   Age is centered within each diagnosis before regression.
#   Therefore, the intercept tests whether the expected Z-score
#   at that disease group's mean age differs from zero.
# ============================================================

suppressPackageStartupMessages({
  library(gamlss)
  library(dplyr)
  library(splines)
})

# ------------------------------------------------------------
# 0) User settings
# ------------------------------------------------------------
BestModelPath <- oef_output("lifespan_model/data/best_model.rds")
HCPath        <- oef_input("trust_lifespan_hc_zscores.csv")
PatientPath   <- oef_input("trust_disease_patients.csv")

OutDir <- oef_output("disease_age_adjustment")
if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)

PhenotypeCol    <- "phenotype"
DiagnosisCol    <- "Diagnosis"
MinSamplesPerDx <- 5
Zthr            <- 1.96

# Same GAMLSS specification used in the original disease script
DF_mu    <- 3
DF_sigma <- 3

# ------------------------------------------------------------
# 1) Load data and basic checks
# ------------------------------------------------------------
# The supplied analysis refits a fixed BCTo model; this file is provenance-only.
best_model <- if (file.exists(BestModelPath)) readRDS(BestModelPath) else NULL
M_HC         <- read.csv(HCPath)
patient_data <- read.csv(PatientPath)

required_patient_cols <- c(PhenotypeCol, DiagnosisCol, "Age", "Sex")
missing_patient_cols <- setdiff(required_patient_cols, names(patient_data))
if (length(missing_patient_cols) > 0) {
  stop("Missing patient columns: ", paste(missing_patient_cols, collapse = ", "))
}

required_hc_cols <- c("phenotype", "Age", "Sex", "SiteID")
missing_hc_cols <- setdiff(required_hc_cols, names(M_HC))
if (length(missing_hc_cols) > 0) {
  stop("Missing HC columns: ", paste(missing_hc_cols, collapse = ", "))
}

cat("\n=== Data overview ===\n")
cat("HC N:", nrow(M_HC), "\n")
cat("Patient N:", nrow(patient_data), "\n")
cat(sprintf("Patient age range: %.2f - %.2f years\n",
            min(patient_data$Age, na.rm = TRUE),
            max(patient_data$Age, na.rm = TRUE)))
cat("Diagnosis counts:\n")
print(table(patient_data[[DiagnosisCol]], useNA = "ifany"))

# ------------------------------------------------------------
# 2) Refit the same clean GAMLSS normative model on HC
# ------------------------------------------------------------
cat("\n=== Refitting clean GAMLSS normative model on HC ===\n")
con <- gamlss.control(n.cyc = 200)

clean_model <- oef_fit_gamlss(
  phenotype ~ bs(Age, df = DF_mu) * Sex + random(as.factor(SiteID)),
  sigma.fo  = ~ bs(Age, df = DF_sigma) + Sex,
  nu.fo     = ~ 1,
  tau.fo    = ~ 1,
  family    = BCTo,
  data      = M_HC,
  control   = con
)

# ------------------------------------------------------------
# ------------------------------------------------------------
newdata_pat <- data.frame(
  Age = patient_data$Age,
  Sex = patient_data$Sex
)

train_site_levels <- levels(factor(M_HC$SiteID))
site_dummy <- train_site_levels[1]
newdata_pat$SiteID <- factor(
  rep(site_dummy, nrow(newdata_pat)),
  levels = train_site_levels
)

cat("\n=== Predicting normative distribution parameters ===\n")
pred_all <- oef_predict_all(random = "zero", clean_model, newdata = newdata_pat)

y <- patient_data[[PhenotypeCol]]
cumulative_prob <- pBCTo(
  y,
  mu    = pred_all$mu,
  sigma = pred_all$sigma,
  nu    = pred_all$nu,
  tau   = pred_all$tau
)

# Avoid +/-Inf after qnorm
cumulative_prob <- pmin(pmax(cumulative_prob, 1e-10), 1 - 1e-10)
patient_data$z_score <- qnorm(cumulative_prob)

write.csv(
  patient_data,
  file.path(OutDir, "disease_zscores_for_age_adjusted_analysis.csv"),
  row.names = FALSE
)

cat(sprintf(
  "Patient Z-score: mean=%.3f, SD=%.3f, range=[%.3f, %.3f]\n",
  mean(patient_data$z_score, na.rm = TRUE),
  sd(patient_data$z_score, na.rm = TRUE),
  min(patient_data$z_score, na.rm = TRUE),
  max(patient_data$z_score, na.rm = TRUE)
))

# ------------------------------------------------------------
# 4) Keep diagnosis groups with adequate sample size
# ------------------------------------------------------------
dx_counts <- table(patient_data[[DiagnosisCol]])
valid_dx <- names(dx_counts[dx_counts >= MinSamplesPerDx])

patient_f <- patient_data %>%
  filter(.data[[DiagnosisCol]] %in% valid_dx) %>%
  mutate(Diagnosis = factor(.data[[DiagnosisCol]], levels = valid_dx))

cat("\n=== Included diagnosis groups ===\n")
print(table(patient_f$Diagnosis))

# ------------------------------------------------------------
# 5) Per-disease age-adjusted regression
#
# Model:
#   z_score = beta0 + beta1 * Age_centered + error
#
# Age_centered = Age - mean(Age within that diagnosis)
#
# Interpretation:
#   beta0 = age-adjusted expected Z-score at the mean age of
#           that disease group.
#   H0: beta0 = 0 tests whether the disease group still deviates
#       from its age/sex-specific normative expectation after
#       explicitly adding chronological age as a covariate.
# ------------------------------------------------------------
results_list <- list()

for (dx in levels(patient_f$Diagnosis)) {

  diag_data <- patient_f %>%
    filter(Diagnosis == dx) %>%
    filter(is.finite(z_score), is.finite(Age))

  N <- nrow(diag_data)
  if (N < MinSamplesPerDx) next

  age_mean <- mean(diag_data$Age)
  age_sd   <- sd(diag_data$Age)
  age_min  <- min(diag_data$Age)
  age_max  <- max(diag_data$Age)

  diag_data <- diag_data %>%
    mutate(Age_centered = Age - age_mean)

  cat("\n----------------------------------------\n")
  cat(sprintf("Diagnosis: %s | N=%d | Mean age=%.2f\n", dx, N, age_mean))

  # Original one-sample t-test for side-by-side comparison
  t_res <- t.test(diag_data$z_score, mu = 0)

  # Reviewer-requested age-adjusted sensitivity model
  fit <- lm(z_score ~ Age_centered, data = diag_data)
  ct  <- coef(summary(fit))
  ci  <- confint(fit, level = 0.95)

  # Intercept: adjusted disease deviation at mean age
  int_est <- ct["(Intercept)", "Estimate"]
  int_se  <- ct["(Intercept)", "Std. Error"]
  int_t   <- ct["(Intercept)", "t value"]
  int_p   <- ct["(Intercept)", "Pr(>|t|)"]
  int_lcl <- ci["(Intercept)", 1]
  int_ucl <- ci["(Intercept)", 2]

  # Residual age association after normative Z-scoring
  age_est <- ct["Age_centered", "Estimate"]
  age_se  <- ct["Age_centered", "Std. Error"]
  age_t   <- ct["Age_centered", "t value"]
  age_p   <- ct["Age_centered", "Pr(>|t|)"]
  age_lcl <- ci["Age_centered", 1]
  age_ucl <- ci["Age_centered", 2]

  results_list[[as.character(dx)]] <- data.frame(
    Diagnosis = as.character(dx),
    N = N,
    Age_mean = age_mean,
    Age_SD = age_sd,
    Age_min = age_min,
    Age_max = age_max,

    Mean_Z_unadjusted = mean(diag_data$z_score),
    SD_Z = sd(diag_data$z_score),
    Ttest_t = unname(t_res$statistic),
    Ttest_df = unname(t_res$parameter),
    Ttest_p = t_res$p.value,

    Age_adjusted_Z = int_est,
    Age_adjusted_SE = int_se,
    Age_adjusted_CI95_L = int_lcl,
    Age_adjusted_CI95_U = int_ucl,
    Age_adjusted_t = int_t,
    Age_adjusted_df = df.residual(fit),
    Age_adjusted_p = int_p,

    Age_beta_per_year = age_est,
    Age_beta_SE = age_se,
    Age_beta_CI95_L = age_lcl,
    Age_beta_CI95_U = age_ucl,
    Age_beta_t = age_t,
    Age_beta_p = age_p,

    Abnormal_percent = mean(abs(diag_data$z_score) > Zthr) * 100,
    stringsAsFactors = FALSE
  )

  cat(sprintf(
    "Original: mean Z=%.3f, p=%.4g | Age-adjusted: Z=%.3f [%.3f, %.3f], p=%.4g | Age beta=%.4f, p=%.4g\n",
    mean(diag_data$z_score), t_res$p.value,
    int_est, int_lcl, int_ucl, int_p,
    age_est, age_p
  ))
}

results_df <- bind_rows(results_list)

# Multiple-testing columns are included for transparency.
# The analysis report can report raw P values if Fig. 2f used raw P values,
# while BH-FDR columns allow direct sensitivity checking across diagnoses.
results_df <- results_df %>%
  mutate(
    Ttest_p_FDR_BH = p.adjust(Ttest_p, method = "BH"),
    Age_adjusted_p_FDR_BH = p.adjust(Age_adjusted_p, method = "BH"),
    Age_beta_p_FDR_BH = p.adjust(Age_beta_p, method = "BH")
  ) %>%
  arrange(desc(Age_adjusted_Z))

# ------------------------------------------------------------
# 6) Save full numeric results table
# ------------------------------------------------------------
full_out <- file.path(OutDir, "disease_age_adjusted_regression_full.csv")
write.csv(results_df, full_out, row.names = FALSE)

# ------------------------------------------------------------
# 7) Save compact publication/reviewer table
# ------------------------------------------------------------
format_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

compact_df <- results_df %>%
  transmute(
    Diagnosis,
    N,
    `Age, mean ± SD` = sprintf("%.1f ± %.1f", Age_mean, Age_SD),
    `Age range` = sprintf("%.1f–%.1f", Age_min, Age_max),
    `Original mean Z` = sprintf("%.3f", Mean_Z_unadjusted),
    `Original t-test P` = format_p(Ttest_p),
    `Age-adjusted Z (95% CI)` = sprintf(
      "%.3f (%.3f, %.3f)",
      Age_adjusted_Z, Age_adjusted_CI95_L, Age_adjusted_CI95_U
    ),
    `Age-adjusted P` = format_p(Age_adjusted_p),
    `Age-adjusted FDR P` = format_p(Age_adjusted_p_FDR_BH),
    `Age beta/year (95% CI)` = sprintf(
      "%.4f (%.4f, %.4f)",
      Age_beta_per_year, Age_beta_CI95_L, Age_beta_CI95_U
    ),
    `Age-effect P` = format_p(Age_beta_p)
  )

compact_out <- file.path(OutDir, "disease_age_adjusted_regression_compact.csv")
write.csv(compact_df, compact_out, row.names = FALSE, fileEncoding = "UTF-8")

# Optional Excel workbook if writexl is already installed.
if (requireNamespace("writexl", quietly = TRUE)) {
  xlsx_out <- file.path(OutDir, "disease_age_adjusted_regression_tables.xlsx")
  writexl::write_xlsx(
    list(
      Full_numeric_results = results_df,
      Compact_reviewer_table = compact_df
    ),
    path = xlsx_out
  )
  cat("Saved Excel workbook to:\n", xlsx_out, "\n")
} else {
  cat("\nPackage 'writexl' not installed; CSV tables were saved normally.\n")
}

cat("\n============================================================\n")
cat("Age-adjusted disease sensitivity analysis completed.\n")
cat("Full results:    ", full_out, "\n")
cat("Compact results: ", compact_out, "\n")
cat("============================================================\n")

print(compact_df)
