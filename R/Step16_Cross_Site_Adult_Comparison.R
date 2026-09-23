# Step16_Cross_Site_Adult_Comparison.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 16
# Scientific and validation caveats: docs/ANALYSIS_NOTES.md
source("R/lib/oef_utils.R")
oef_start_step("Step16")

# ============================================================
# Site 1 vs traveling-validation Sites 8 and 12
# Age-overlapping adult comparison: 23–33 years
# ============================================================

library(readxl)

# ------------------------------------------------------------
# 1. Read FINAL HC dataset
# ------------------------------------------------------------

file <- oef_input("trust_hc_hct_ya_sensitivity.xlsx")

OutDir <- file.path(oef_output(), "cross_site_adults")
dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)
dat_all <- read_excel(file, sheet = "Sheet1")
if (!("Sex" %in% names(dat_all)) && "Sex(M0F1)" %in% names(dat_all)) {
  dat_all$Sex <- dat_all[["Sex(M0F1)"]]
}
oef_require_columns(dat_all, c("SiteID", "Age", "Sex", "T2", "OEF"))

# Check column names
print(names(dat_all))


# ------------------------------------------------------------
# 2. Select Sites 1, 8, and 12
#    Restrict to common adult age range: 23–33 years
# ------------------------------------------------------------

dat <- dat_all[
  dat_all$SiteID %in% c(1, 8, 12) &
    dat_all$Age >= 23 &
    dat_all$Age <= 33,
]

# Sex retains the input coding, 0=Male and 1=Female.

# Site indicator:
#   Site1 = 1: Site 1
#   Site1 = 0: Sites 8 and 12
dat$Site1 <- ifelse(dat$SiteID == 1, 1, 0)

# Calculate R2 from venous blood T2
# T2 is stored in milliseconds
# R2 (s^-1) = 1/T2(s) = 1000/T2(ms)
if (any(!is.na(dat$T2) & (!is.finite(dat$T2) | dat$T2 <= 0))) stop("T2 must be positive milliseconds.")
dat$R2_per_s <- 1000 / dat$T2


# ------------------------------------------------------------
# 3. Check sample composition
# ------------------------------------------------------------

cat("\nSample size by site:\n")
print(table(dat$SiteID))

cat("\nAge distribution by site:\n")
print(
  aggregate(
    Age ~ SiteID,
    data = dat,
    FUN = function(x)
      c(
        N = length(x),
        mean = mean(x),
        SD = sd(x),
        min = min(x),
        max = max(x)
      )
  )
)

cat("\nSex distribution by site:\n")
print(table(dat$SiteID, dat$Sex))


# ------------------------------------------------------------
# 4. OEF comparison
#
# Model:
# OEF ~ Site1 + Age + Sex
# ------------------------------------------------------------

model_oef <- lm(
  OEF ~ Site1 + Age + Sex,
  data = dat
)

cat("\n====================================\n")
cat("OEF MODEL\n")
cat("====================================\n")

print(summary(model_oef))

# Coefficient
beta_oef <- coef(model_oef)["Site1"]

# 95% CI
ci_oef <- confint(model_oef, "Site1", level = 0.95)

# P value
p_oef <- summary(model_oef)$coefficients["Site1", "Pr(>|t|)"]

cat("\nAdjusted OEF difference: Site 1 - Sites 8/12\n")

cat(
  sprintf(
    "Beta = %.6f (OEF fraction)\n",
    beta_oef
  )
)

cat(
  sprintf(
    "Difference = %.3f percentage points\n",
    beta_oef * 100
  )
)

cat(
  sprintf(
    "95%% CI = [%.3f, %.3f] percentage points\n",
    ci_oef[1] * 100,
    ci_oef[2] * 100
  )
)

cat(
  sprintf(
    "P = %.4f\n",
    p_oef
  )
)


# ------------------------------------------------------------
# 5. R2 comparison
#
# Model:
# R2 ~ Site1 + Age + Sex
# ------------------------------------------------------------

model_r2 <- lm(
  R2_per_s ~ Site1 + Age + Sex,
  data = dat
)

cat("\n====================================\n")
cat("R2 MODEL\n")
cat("====================================\n")

print(summary(model_r2))

# Coefficient
beta_r2 <- coef(model_r2)["Site1"]

# 95% CI
ci_r2 <- confint(model_r2, "Site1", level = 0.95)

# P value
p_r2 <- summary(model_r2)$coefficients["Site1", "Pr(>|t|)"]

cat("\nAdjusted R2 difference: Site 1 - Sites 8/12\n")

cat(
  sprintf(
    "Difference = %.3f s^-1\n",
    beta_r2
  )
)

cat(
  sprintf(
    "95%% CI = [%.3f, %.3f] s^-1\n",
    ci_r2[1],
    ci_r2[2]
  )
)

cat(
  sprintf(
    "P = %.4f\n",
    p_r2
  )
)

# Persist the exact fitted contrast, not only printed rounded values.
write.csv(data.frame(
  Outcome = c("OEF_fraction", "R2_per_s"),
  N = c(nobs(model_oef), nobs(model_r2)),
  Estimate = c(unname(beta_oef), unname(beta_r2)),
  CI_low = c(ci_oef[1], ci_r2[1]), CI_high = c(ci_oef[2], ci_r2[2]),
  P = c(unname(p_oef), unname(p_r2))
), file.path(OutDir, "site1_vs_sites8_12_adjusted_contrasts.csv"), row.names = FALSE)
write.csv(as.data.frame(table(dat$SiteID)), file.path(OutDir, "site_counts.csv"), row.names = FALSE)
