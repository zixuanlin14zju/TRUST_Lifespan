# Step8_Hct_Sensitivity_Normative.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 8
source("R/lib/oef_utils.R")
oef_start_step("Step8")

# ============================================================
# TRUST lifespan OEF hematocrit sensitivity: normative trajectories
# Reference OEF plus RegionMatched, Zierk, MahlknechtMean and CohortEmpirical.
# The reference model selects df(mu) by BIC; all four sensitivity models
# use that same df. The distribution, sigma model, and prediction method
# are unchanged from the existing lifespan analysis.
#
# Input: trust_hc_hct_ya_sensitivity.xlsx, Sheet1.
# Outputs: model_selection.csv, sex_interaction_LRT.csv,
# trajectory_by_age.csv, growth_rate_by_age.csv,
# sensitivity_vs_original.csv, change_age20_to_70.csv, fitted models.
# Overall/female/male figures are generated from these outputs by Step19.
# ============================================================

# ============================================================
# 00) Packages
# ============================================================

required_packages <- c(
  "gamlss",
  "pracma",
  "dplyr",
  "splines",
  "readxl"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install the following package(s) first: ",
    paste(missing_packages, collapse = ", "),
    "\nExample: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(gamlss)
  library(pracma)
  library(dplyr)
  library(splines)
  library(readxl)
})

# ============================================================
# 01) User settings
# ============================================================

HCFile <- oef_input("trust_hc_hct_ya_sensitivity.xlsx")

OutDir <- file.path(oef_output(), "GAMLSS_HctSensitivity")
dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)

# Main-model specification
CandidateDfMu <- 3:5
DfSigma       <- 3

# Age grid for trajectory output
AgeStep <- 0.1

# ============================================================
# 02) Scenario definitions
# ============================================================

# "Original" is the submitted/reference OEF and is not counted as a
# sensitivity scenario. The four requested sensitivity analyses follow it.
ScenarioMap <- c(
  Original        = "OEF",
  RegionMatched   = "OEF_RegionMatched",
  Zierk           = "OEF_Zierk",
  MahlknechtMean  = "OEF_MahlMean",
  CohortEmpirical = "OEF_CohortEmpirical"
)

SensitivityScenarios <- setdiff(names(ScenarioMap), "Original")

# ============================================================
# 03) Utility functions
# ============================================================

normalize_sex <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    out <- suppressWarnings(as.numeric(x))
    if (!all(na.omit(unique(out)) %in% c(0, 1))) {
      stop("Numeric Sex values are not coded as 0/1.")
    }
    return(as.integer(out))
  }

  s <- tolower(trimws(as.character(x)))

  out <- ifelse(
    s %in% c("1", "f", "female", "woman", "w"),
    1L,
    ifelse(
      s %in% c("0", "m", "male", "man"),
      0L,
      NA_integer_
    )
  )

  if (any(is.na(out) & !is.na(x))) {
    stop("Could not map all Sex values to 0=Male / 1=Female.")
  }

  out
}

prepare_input <- function(path, dataset_label) {

  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }

  dat <- as.data.frame(
    readxl::read_excel(path, sheet = "Sheet1"),
    stringsAsFactors = FALSE
  )

  # Current Excel uses Sex(M0F1); retain support for an existing Sex column.
  if ("Sex(M0F1)" %in% names(dat)) {
    dat$Sex <- normalize_sex(dat[["Sex(M0F1)"]])
  } else if ("Sex" %in% names(dat)) {
    dat$Sex <- normalize_sex(dat[["Sex"]])
  } else {
    stop(dataset_label, ": no Sex or Sex(M0F1) column found.")
  }

  required <- c("Age", "SiteID", unname(ScenarioMap))
  missing <- setdiff(required, names(dat))

  if (length(missing) > 0) {
    stop(
      dataset_label,
      ": missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }

  dat$Age <- suppressWarnings(as.numeric(dat$Age))
  dat$SiteID <- as.character(dat$SiteID)

  # Use QC=PASS when that column exists.
  if ("QC" %in% names(dat)) {
    n_before <- nrow(dat)
    keep_qc <- is.na(dat$QC) | toupper(trimws(as.character(dat$QC))) == "PASS"
    dat <- dat[keep_qc, , drop = FALSE]

    cat(
      sprintf(
        "%s: QC filter retained %d/%d rows.\n",
        dataset_label,
        nrow(dat),
        n_before
      )
    )
  }

  dat
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok])
}

safe_rmse <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 1) return(NA_real_)
  sqrt(mean((x[ok] - y[ok])^2))
}

# ============================================================
# 04) Load HC data
# ============================================================
HC <- prepare_input(HCFile, "HC")
available <- sapply(unname(ScenarioMap), function(nm) is.finite(HC[[nm]]))
coverage <- data.frame(
  Cohort = "HC", Scenario = names(ScenarioMap),
  N_available = colSums(available),
  N_complete_all_scenarios = sum(rowSums(available) == ncol(available))
)
write.csv(coverage, file.path(OutDir, "scenario_input_coverage.csv"), row.names = FALSE)
if (any(coverage$N_available != coverage$N_complete_all_scenarios)) {
  warning("Hct scenarios have unequal available rows. Scenario-specific filtering is retained; inspect scenario_input_coverage.csv.")
}
cat("\nHct sensitivity: HC rows after QC =", nrow(HC), "\n")
cat("Scenarios:", paste(names(ScenarioMap), collapse = ", "), "\n")

# Common age range, based on HC.
CommonAgeMin <- min(HC$Age, na.rm = TRUE)
CommonAgeMax <- max(HC$Age, na.rm = TRUE)

AgeGrid <- seq(
  CommonAgeMin,
  CommonAgeMax,
  by = AgeStep
)

# ============================================================
# 05) GAMLSS model helpers
# ============================================================

gamlss_con <- gamlss.control(n.cyc = 200, trace = FALSE)

fit_one_gamlss <- function(dat, df_mu, include_interaction = TRUE) {

  if (include_interaction) {
    mu_formula <- phenotype ~
      bs(Age, df = df_mu) * Sex +
      random(as.factor(SiteID))
  } else {
    mu_formula <- phenotype ~
      bs(Age, df = df_mu) +
      Sex +
      random(as.factor(SiteID))
  }

  oef_fit_gamlss(
    formula  = mu_formula,
    sigma.fo = ~ bs(Age, df = DfSigma) + Sex,
    nu.fo    = ~ 1,
    tau.fo   = ~ 1,
    family   = BCTo,
    data     = dat,
    control  = gamlss_con
  )
}

select_df_by_bic <- function(dat) {

  models <- vector("list", length(CandidateDfMu))
  names(models) <- as.character(CandidateDfMu)

  audit <- data.frame(
    df_mu = CandidateDfMu,
    converged = FALSE,
    BIC = NA_real_
  )

  for (i in seq_along(CandidateDfMu)) {

    d <- CandidateDfMu[i]

    cat("  Candidate df(mu) =", d, "...\n")

    fit <- tryCatch(
      fit_one_gamlss(dat, d, include_interaction = TRUE),
      error = function(e) {
        message("    Fit failed: ", e$message)
        NULL
      }
    )

    models[i] <- list(fit)

    if (!is.null(fit)) {
      audit$converged[i] <- isTRUE(fit$converged)
      audit$BIC[i] <- fit$sbc
    }
  }

  valid <- which(
    audit$converged &
    is.finite(audit$BIC)
  )

  if (length(valid) == 0) {
    stop("No candidate GAMLSS model converged.")
  }

  best_i <- valid[which.min(audit$BIC[valid])]

  list(
    best_df = audit$df_mu[best_i],
    best_model = models[[best_i]],
    audit = audit
  )
}

fit_locked_df <- function(dat, df_mu) {

  fit <- fit_one_gamlss(
    dat,
    df_mu,
    include_interaction = TRUE
  )

  if (!isTRUE(fit$converged)) {
    warning(
      "Locked-df model did not report convergence at df(mu)=",
      df_mu
    )
  }

  fit
}

sex_interaction_lrt <- function(dat, full_model, df_mu) {

  reduced_model <- tryCatch(
    fit_one_gamlss(
      dat,
      df_mu,
      include_interaction = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(reduced_model)) {
    return(
      data.frame(
        LR = NA_real_,
        df_diff = NA_real_,
        P = NA_real_,
        reduced_converged = FALSE
      )
    )
  }

  ll_full <- logLik(full_model)
  ll_red  <- logLik(reduced_model)

  LR <- 2 * (as.numeric(ll_full) - as.numeric(ll_red))
  df_diff <- attr(ll_full, "df") - attr(ll_red, "df")

  P <- if (
    is.finite(LR) &&
    is.finite(df_diff) &&
    df_diff > 0
  ) {
    pchisq(max(LR, 0), df = df_diff, lower.tail = FALSE)
  } else {
    NA_real_
  }

  data.frame(
    LR = LR,
    df_diff = df_diff,
    P = P,
    reduced_converged = isTRUE(reduced_model$converged)
  )
}

# ============================================================
# 06) Prediction helpers
# ============================================================

predict_parameters <- function(model, age, sex, site_levels) {

  nd <- data.frame(
    Age = age,
    Sex = sex,
    SiteID = factor(
      rep(site_levels[1], length(age)),
      levels = site_levels
    )
  )

  oef_predict_all(random = "zero", 
    model,
    newdata = nd
  )
}

predict_trajectory <- function(model, scenario, age_grid, site_levels) {

  pa_f <- predict_parameters(
    model,
    age_grid,
    sex = 1,
    site_levels = site_levels
  )

  pa_m <- predict_parameters(
    model,
    age_grid,
    sex = 0,
    site_levels = site_levels
  )

  # Median centile from fitted BCTo distribution.
  med_f <- qBCTo(
    0.5,
    mu = pa_f$mu,
    sigma = pa_f$sigma,
    nu = pa_f$nu,
    tau = pa_f$tau
  )

  med_m <- qBCTo(
    0.5,
    mu = pa_m$mu,
    sigma = pa_m$sigma,
    nu = pa_m$nu,
    tau = pa_m$tau
  )

  med_all <- (med_f + med_m) / 2

  # Keep the original submitted pipeline's growth-rate definition:
  mu_all <- (pa_f$mu + pa_m$mu) / 2

  growth_f <- pracma::gradient(pa_f$mu, age_grid)
  growth_m <- pracma::gradient(pa_m$mu, age_grid)
  growth_all <- pracma::gradient(mu_all, age_grid)

  trajectory <- data.frame(
    Scenario = scenario,
    Age = age_grid,
    Median_Overall = med_all,
    Median_Female = med_f,
    Median_Male = med_m
  )

  growth <- data.frame(
    Scenario = scenario,
    Age = age_grid,
    Growth_Overall = growth_all,
    Growth_Female = growth_f,
    Growth_Male = growth_m
  )

  list(
    trajectory = trajectory,
    growth = growth
  )
}

# ============================================================
# 07) Fit reference model first and determine locked df
# ============================================================

ModelObjects <- list()
TrajectoryList <- list()
GrowthList <- list()
ModelSelectionList <- list()
SexLRTList <- list()

ReferenceDf <- NA_integer_

for (scenario in names(ScenarioMap)) {

  phenotype_col <- ScenarioMap[[scenario]]

  cat("\n============================================================\n")
  cat("SCENARIO:", scenario, "\n")
  cat("Phenotype column:", phenotype_col, "\n")
  cat("============================================================\n")

  dat <- HC

  dat$phenotype <- suppressWarnings(
    as.numeric(dat[[phenotype_col]])
  )

  # Keep ONLY variables actually used by the GAMLSS model.
  dat <- dat[, c("Age", "Sex", "SiteID", "phenotype"), drop = FALSE]

  dat <- dat[
    is.finite(dat$Age) &
    !is.na(dat$Sex) &
    !is.na(dat$SiteID) &
    is.finite(dat$phenotype),
    ,
    drop = FALSE
  ]

  dat <- na.omit(dat)
  dat$SiteID <- as.factor(dat$SiteID)

  cat("N used:", nrow(dat), "\n")
  cat("Sites:", nlevels(dat$SiteID), "\n")

  if (scenario == "Original") {
    selection <- select_df_by_bic(dat)
    ReferenceDf <- selection$best_df
    model <- selection$best_model
    scenario_audit <- selection$audit
    scenario_audit$SelectionRole <- "BIC_selection_on_Original"
    cat("  Reference df(mu):", ReferenceDf, "\n")
  } else {
    if (!is.finite(ReferenceDf)) stop("Fit Original before sensitivity scenarios.")
    model <- fit_locked_df(dat, ReferenceDf)
    scenario_audit <- data.frame(
      df_mu = ReferenceDf,
      converged = isTRUE(model$converged),
      BIC = model$sbc,
      SelectionRole = "df_locked_from_Original"
    )
  }
  scenario_audit$Scenario <- scenario
  scenario_audit$N <- nrow(dat)
  ModelSelectionList[[scenario]] <- scenario_audit

  ModelObjects[[scenario]] <- model

  saveRDS(
    model,
    file.path(
      OutDir,
      paste0("best_model_", scenario, ".rds")
    )
  )

  fit_df_used <- ReferenceDf

  # Sex interaction test, using the same df as the final scenario model.
  sex_lrt <- sex_interaction_lrt(
    dat,
    model,
    fit_df_used
  )

  sex_lrt$Scenario <- scenario
  sex_lrt$df_mu_used <- fit_df_used
  sex_lrt$N <- nrow(dat)

  SexLRTList[[scenario]] <- sex_lrt

  cat(
    "  Sex trajectory interaction LRT P =",
    format.pval(sex_lrt$P, digits = 4),
    "\n"
  )

  site_levels <- levels(dat$SiteID)

  pr <- predict_trajectory(
    model,
    scenario,
    AgeGrid,
    site_levels
  )

  TrajectoryList[[scenario]] <- pr$trajectory
  GrowthList[[scenario]] <- pr$growth
}

# ============================================================
# 08) Consolidate model-selection and trajectory outputs
# ============================================================

ModelSelection <- bind_rows(ModelSelectionList) %>%
  select(
    Scenario,
    N,
    df_mu,
    converged,
    BIC,
    SelectionRole
  ) %>%
  arrange(
    factor(Scenario, levels = names(ScenarioMap)),
    df_mu
  )

write.csv(
  ModelSelection,
  file.path(OutDir, "model_selection.csv"),
  row.names = FALSE
)

SexLRT <- bind_rows(SexLRTList) %>%
  select(
    Scenario,
    N,
    df_mu_used,
    LR,
    df_diff,
    P,
    reduced_converged
  )

write.csv(
  SexLRT,
  file.path(OutDir, "sex_interaction_LRT.csv"),
  row.names = FALSE
)

Trajectory <- bind_rows(TrajectoryList)

Growth <- bind_rows(GrowthList)

write.csv(
  Trajectory,
  file.path(OutDir, "trajectory_by_age.csv"),
  row.names = FALSE
)

write.csv(
  Growth,
  file.path(OutDir, "growth_rate_by_age.csv"),
  row.names = FALSE
)

# ============================================================
# 09) Quantitative comparison with Original trajectory
# ============================================================

ref_traj <- Trajectory %>%
  filter(Scenario == "Original") %>%
  arrange(Age)

ref_growth <- Growth %>%
  filter(Scenario == "Original") %>%
  arrange(Age)

comparison_rows <- list()

for (scenario in SensitivityScenarios) {

  tr <- Trajectory %>%
    filter(Scenario == scenario) %>%
    arrange(Age)

  gr <- Growth %>%
    filter(Scenario == scenario) %>%
    arrange(Age)

  # Values at clinically/interpretable ages.
  nearest_value <- function(age_target, x, y) {
    idx <- which.min(abs(x - age_target))
    y[idx]
  }

  med20 <- nearest_value(
    20,
    tr$Age,
    tr$Median_Overall
  )

  med70 <- nearest_value(
    70,
    tr$Age,
    tr$Median_Overall
  )

  ref20 <- nearest_value(
    20,
    ref_traj$Age,
    ref_traj$Median_Overall
  )

  ref70 <- nearest_value(
    70,
    ref_traj$Age,
    ref_traj$Median_Overall
  )

  comparison_rows[[scenario]] <- data.frame(
    Scenario = scenario,

    Trajectory_r_vs_Original =
      safe_cor(
        tr$Median_Overall,
        ref_traj$Median_Overall
      ),

    Trajectory_RMSE_vs_Original =
      safe_rmse(
        tr$Median_Overall,
        ref_traj$Median_Overall
      ),

    Growth_r_vs_Original =
      safe_cor(
        gr$Growth_Overall,
        ref_growth$Growth_Overall
      ),

    Growth_RMSE_vs_Original =
      safe_rmse(gr$Growth_Overall, ref_growth$Growth_Overall),
    OEF20 = med20,
    OEF70 = med70,
    Delta20to70 = med70 - med20,
    Original_Delta20to70 = ref70 - ref20,
    PercentChangeVsOriginal =
      if (is.finite(ref70 - ref20) && ref70 != ref20) {
        100 * ((med70 - med20) - (ref70 - ref20)) / abs(ref70 - ref20)
      } else NA_real_
  )
}

SensitivitySummary <- bind_rows(comparison_rows)

SensitivitySummary$Trajectory_one_minus_r <- 1 - SensitivitySummary$Trajectory_r_vs_Original
SensitivitySummary$Growth_one_minus_r <- 1 - SensitivitySummary$Growth_r_vs_Original

write.csv(
  SensitivitySummary,
  file.path(OutDir, "sensitivity_vs_original.csv"),
  row.names = FALSE
)

# ============================================================
# 10) Age 20-70 summary, including the reference model
# ============================================================
change_rows <- lapply(names(ScenarioMap), function(scenario) {
  tr <- Trajectory[Trajectory$Scenario == scenario, ]
  if (min(tr$Age) > 20 || max(tr$Age) < 70) {
    stop("The fitted age grid does not cover both 20 and 70 years.")
  }
  i20 <- which.min(abs(tr$Age - 20))
  i70 <- which.min(abs(tr$Age - 70))
  data.frame(
    Scenario = scenario,
    OEF_age20 = tr$Median_Overall[i20],
    OEF_age70 = tr$Median_Overall[i70],
    Delta20to70 = tr$Median_Overall[i70] - tr$Median_Overall[i20]
  )
})
ChangeSummary <- bind_rows(change_rows)
ref_delta <- ChangeSummary$Delta20to70[ChangeSummary$Scenario == "Original"]
ChangeSummary$PercentChangeVsOriginal <- if (is.finite(ref_delta) && ref_delta != 0) {
  100 * (ChangeSummary$Delta20to70 - ref_delta) / abs(ref_delta)
} else NA_real_
write.csv(ChangeSummary, file.path(OutDir, "change_age20_to_70.csv"), row.names = FALSE)

cat("\nHC Hct sensitivity complete. Reference df(mu):", ReferenceDf, "\n")
print(SexLRT)
print(SensitivitySummary)
cat("Outputs:", OutDir, "\nRun Step19 after Step9 for the overall/female/male figures.\n")
