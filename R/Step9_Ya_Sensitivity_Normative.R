# Step9_Ya_Sensitivity_Normative.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 9
source("R/lib/oef_utils.R")
oef_start_step("Step9")

# ============================================================
# TRUST lifespan OEF arterial-oxygenation (Ya) sensitivity
# GAMLSS BCTo normative modeling -- FINAL HYBRID PEDIATRIC/ADULT DESIGN
#
# Primary/reference OEF:
#   Original submitted OEF. Site 1 neonates used measured Ya;
#   most other pediatric participants and adults used Ya = 0.98.
#
# Sensitivity 1: LuHybrid
#   Pediatric (<18 y): use individual measured Ya whenever available;
#                      if unavailable, extrapolate the Lu age equation.
#   Adult (>=18 y):    use the Lu age equation for all participants.
#   Lu et al. 2011: Ya(%) = 99.06 - 0.02 * Age(years)
#
# Sensitivity 2: CrapoHybrid
#   Pediatric (<18 y): use individual measured Ya whenever available;
#                      if unavailable, extrapolate the Crapo age equation.
#   Adult (>=18 y):    use the Crapo age equation for all participants.
#   Crapo et al. 1999: SaO2(%) = 97.66 - 0.0296 * Age(years)
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
    paste(missing_packages, collapse = ", ")
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

OutDir <- file.path(oef_output(), "GAMLSS_YaSensitivity_Hybrid_FINAL")
dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)

CandidateDfMu <- 3:5
DfSigma <- 3
AgeStep <- 0.1

# ============================================================
# 02) Utility functions
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

prepare_hc <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)

  dat <- as.data.frame(
    readxl::read_excel(path, sheet = "Sheet1"),
    stringsAsFactors = FALSE
  )

  if ("Sex(M0F1)" %in% names(dat)) {
    dat$Sex <- normalize_sex(dat[["Sex(M0F1)"]])
  } else if ("Sex" %in% names(dat)) {
    dat$Sex <- normalize_sex(dat[["Sex"]])
  } else {
    stop("No Sex or Sex(M0F1) column found.")
  }

  required <- c(
    "Age",
    "SiteID",
    "OEF",
    "OEF_Ya_LuHybrid",
    "OEF_Ya_CrapoHybrid"
  )

  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop("Missing required column(s): ", paste(missing, collapse = ", "))
  }

  dat$Age <- suppressWarnings(as.numeric(dat$Age))
  dat$SiteID <- as.character(dat$SiteID)

  if ("QC" %in% names(dat)) {
    keep <- is.na(dat$QC) |
      toupper(trimws(as.character(dat$QC))) == "PASS"
    cat("QC filter retained", sum(keep), "/", nrow(dat), "rows.\n")
    dat <- dat[keep, , drop = FALSE]
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
  if (!any(ok)) return(NA_real_)
  sqrt(mean((x[ok] - y[ok])^2))
}

# ============================================================
# 03) GAMLSS helpers
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
    formula = mu_formula,
    sigma.fo = ~ bs(Age, df = DfSigma) + Sex,
    nu.fo = ~ 1,
    tau.fo = ~ 1,
    family = BCTo,
    data = dat,
    control = gamlss_con
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
    cat("    Candidate df(mu) =", d, "...\n")

    fit <- tryCatch(
      fit_one_gamlss(dat, d, include_interaction = TRUE),
      error = function(e) {
        message("      Fit failed: ", e$message)
        NULL
      }
    )

    models[i] <- list(fit)

    if (!is.null(fit)) {
      audit$converged[i] <- isTRUE(fit$converged)
      audit$BIC[i] <- fit$sbc
    }
  }

  valid <- which(audit$converged & is.finite(audit$BIC))
  if (length(valid) == 0) stop("No candidate GAMLSS model converged.")

  best_i <- valid[which.min(audit$BIC[valid])]

  list(
    best_df = audit$df_mu[best_i],
    best_model = models[[best_i]],
    audit = audit
  )
}

sex_interaction_lrt <- function(dat, full_model, df_mu) {
  reduced_model <- tryCatch(
    fit_one_gamlss(dat, df_mu, include_interaction = FALSE),
    error = function(e) NULL
  )

  if (is.null(reduced_model)) {
    return(data.frame(LR = NA_real_, df_diff = NA_real_, P = NA_real_))
  }

  ll_full <- logLik(full_model)
  ll_red <- logLik(reduced_model)

  LR <- 2 * (as.numeric(ll_full) - as.numeric(ll_red))
  df_diff <- attr(ll_full, "df") - attr(ll_red, "df")

  P <- if (is.finite(LR) && is.finite(df_diff) && df_diff > 0) {
    pchisq(max(LR, 0), df = df_diff, lower.tail = FALSE)
  } else {
    NA_real_
  }

  data.frame(LR = LR, df_diff = df_diff, P = P)
}

predict_parameters <- function(model, age, sex, site_levels) {
  nd <- data.frame(
    Age = age,
    Sex = sex,
    SiteID = factor(
      rep(site_levels[1], length(age)),
      levels = site_levels
    )
  )

  oef_predict_all(random = "zero", model, newdata = nd)
}

predict_trajectory <- function(model, scenario, age_grid, site_levels) {
  pf <- predict_parameters(model, age_grid, 1, site_levels)
  pm <- predict_parameters(model, age_grid, 0, site_levels)

  med_f <- qBCTo(
    0.5,
    mu = pf$mu,
    sigma = pf$sigma,
    nu = pf$nu,
    tau = pf$tau
  )

  med_m <- qBCTo(
    0.5,
    mu = pm$mu,
    sigma = pm$sigma,
    nu = pm$nu,
    tau = pm$tau
  )

  med_all <- (med_f + med_m) / 2
  mu_all <- (pf$mu + pm$mu) / 2

  data.frame(
    Scenario = scenario,
    Age = age_grid,
    Median_Overall = med_all,
    Median_Female = med_f,
    Median_Male = med_m,
    Growth_Overall = pracma::gradient(mu_all, age_grid),
    Growth_Female = pracma::gradient(pf$mu, age_grid),
    Growth_Male = pracma::gradient(pm$mu, age_grid)
  )
}

# ============================================================
# 04) One analysis-set runner
# ============================================================
run_analysis_set <- function(
    HC,
    analysis_name,
    scenario_map,
    age_filter) {

  this_out <- file.path(OutDir, analysis_name)
  dir.create(this_out, recursive = TRUE, showWarnings = FALSE)

  cat("\n============================================================\n")
  cat("ANALYSIS SET:", analysis_name, "\n")
  cat("============================================================\n")

  base <- HC[age_filter(HC$Age), , drop = FALSE]

  # Require complete values for ALL scenarios so every model in this
  # analysis set uses exactly the same participants.
  required_scenario_cols <- unname(scenario_map)
  complete_all <- is.finite(base$Age) &
    !is.na(base$Sex) &
    !is.na(base$SiteID)

  for (col in required_scenario_cols) {
    complete_all <- complete_all & is.finite(suppressWarnings(as.numeric(base[[col]])))
  }

  base <- base[complete_all, , drop = FALSE]

  cat("Common N across scenarios:", nrow(base), "\n")
  cat("Age range:", min(base$Age), "to", max(base$Age), "\n")
  cat("Sites:", length(unique(base$SiteID)), "\n")

  age_grid <- seq(min(base$Age), max(base$Age), by = AgeStep)

  models <- list()
  model_selection <- list()
  sex_lrt <- list()
  trajectories <- list()

  reference_df <- NA_integer_

  for (scenario in names(scenario_map)) {
    phenotype_col <- scenario_map[[scenario]]

    cat("\n  SCENARIO:", scenario, "\n")
    cat("  Column:", phenotype_col, "\n")

    dat <- base
    dat$phenotype <- suppressWarnings(as.numeric(dat[[phenotype_col]]))

    # gamlss checks the entire data.frame for NA values, so pass only
    # variables actually used by the model.
    dat <- dat[, c("Age", "Sex", "SiteID", "phenotype"), drop = FALSE]
    dat <- na.omit(dat)
    dat$SiteID <- factor(dat$SiteID)

    if (scenario == "Original") {
      selection <- select_df_by_bic(dat)
      reference_df <- selection$best_df
      model <- selection$best_model

      audit <- selection$audit
      audit$SelectionRole <- "BIC_selection_on_Original"

      cat("    Reference BIC-selected df(mu) =", reference_df, "\n")
    } else {
      if (!is.finite(reference_df)) {
        stop("Original scenario must be fitted before sensitivity scenarios.")
      }

      cat("    Locked df(mu) =", reference_df, "\n")
      model <- fit_one_gamlss(
        dat,
        df_mu = reference_df,
        include_interaction = TRUE
      )

      audit <- data.frame(
        df_mu = reference_df,
        converged = isTRUE(model$converged),
        BIC = model$sbc,
        SelectionRole = "df_locked_from_Original"
      )
    }

    audit$Analysis <- analysis_name
    audit$Scenario <- scenario
    audit$N <- nrow(dat)
    model_selection[[scenario]] <- audit

    models[[scenario]] <- model
    saveRDS(model, file.path(this_out, paste0("model_", scenario, ".rds")))

    sx <- sex_interaction_lrt(dat, model, reference_df)
    sx$Analysis <- analysis_name
    sx$Scenario <- scenario
    sx$df_mu_used <- reference_df
    sx$N <- nrow(dat)
    sex_lrt[[scenario]] <- sx

    site_levels <- levels(dat$SiteID)
    trajectories[[scenario]] <- predict_trajectory(
      model,
      scenario,
      age_grid,
      site_levels
    )

    cat("    Sex interaction LRT P =", format.pval(sx$P, digits = 4), "\n")
  }

  model_selection_df <- bind_rows(model_selection)
  sex_lrt_df <- bind_rows(sex_lrt)
  trajectory_df <- bind_rows(trajectories)

  write.csv(
    model_selection_df,
    file.path(this_out, "model_selection.csv"),
    row.names = FALSE
  )

  write.csv(
    sex_lrt_df,
    file.path(this_out, "sex_interaction_LRT.csv"),
    row.names = FALSE
  )

  write.csv(
    trajectory_df,
    file.path(this_out, "trajectory_and_growth_by_age.csv"),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # Scenario concordance vs Original
  # ----------------------------------------------------------
  ref <- trajectory_df[trajectory_df$Scenario == "Original", ]
  comp_list <- list()

  for (scenario in setdiff(names(scenario_map), "Original")) {
    xx <- trajectory_df[trajectory_df$Scenario == scenario, ]

    comp_list[[scenario]] <- data.frame(
      Analysis = analysis_name,
      Scenario = scenario,
      Trajectory_r = safe_cor(ref$Median_Overall, xx$Median_Overall),
      Trajectory_RMSE = safe_rmse(ref$Median_Overall, xx$Median_Overall),
      Growth_r = safe_cor(ref$Growth_Overall, xx$Growth_Overall),
      Growth_RMSE = safe_rmse(ref$Growth_Overall, xx$Growth_Overall)
    )
  }

  concordance_df <- bind_rows(comp_list)
  concordance_df$Trajectory_one_minus_r <- 1 - concordance_df$Trajectory_r
  concordance_df$Growth_one_minus_r <- 1 - concordance_df$Growth_r
  write.csv(
    concordance_df,
    file.path(this_out, "sensitivity_vs_original.csv"),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # 20 -> 70 year change when available
  # ----------------------------------------------------------
  change_list <- list()

  for (scenario in names(scenario_map)) {
    tr <- trajectory_df[trajectory_df$Scenario == scenario, ]

    if (min(tr$Age) <= 20 && max(tr$Age) >= 70) {
      i20 <- which.min(abs(tr$Age - 20))
      i70 <- which.min(abs(tr$Age - 70))

      change_list[[scenario]] <- data.frame(
        Analysis = analysis_name,
        Scenario = scenario,
        OEF_age20 = tr$Median_Overall[i20],
        OEF_age70 = tr$Median_Overall[i70],
        Delta20to70 = tr$Median_Overall[i70] - tr$Median_Overall[i20]
      )
    }
  }

  change_df <- bind_rows(change_list)

  if (nrow(change_df) > 0) {
    ref_delta <- change_df$Delta20to70[change_df$Scenario == "Original"]
    change_df$PercentChangeVsOriginal <-
      100 * (change_df$Delta20to70 - ref_delta) / abs(ref_delta)

    write.csv(
      change_df,
      file.path(this_out, "change_age20_to_70.csv"),
      row.names = FALSE
    )
  }

  # ----------------------------------------------------------
  invisible(list(
    models = models,
    trajectory = trajectory_df,
    concordance = concordance_df,
    changes = change_df
  ))
}

# ============================================================
# 05) Load data
# ============================================================
HC <- prepare_hc(HCFile)

cat("\n============================================================\n")
cat("Ya sensitivity analysis\n")
cat("HC rows after QC:", nrow(HC), "\n")
cat("============================================================\n")

# ============================================================
# 06) Full-lifespan hybrid Ya sensitivity
# ============================================================
HybridMap <- c(
  Original = "OEF",
  LuHybrid = "OEF_Ya_LuHybrid",
  CrapoHybrid = "OEF_Ya_CrapoHybrid"
)

res_hybrid <- run_analysis_set(
  HC = HC,
  analysis_name = "01_FullLifespan_Original_vs_LuHybrid_vs_CrapoHybrid",
  scenario_map = HybridMap,
  age_filter = function(age) is.finite(age)
)

# ============================================================
# 07) Save analysis notes
# ============================================================
notes <- c(
  "TRUST OEF Ya sensitivity analysis -- hybrid pediatric/adult design",
  "",
  "Primary/reference:",
  "  Original submitted OEF. Site 1 neonates used measured Ya; most other pediatric participants and adults used Ya=0.98.",
  "",
  "LuHybrid:",
  "  Pediatric (<18 y): measured individual Ya when available; otherwise Lu equation extrapolation.",
  "  Adult (>=18 y): Lu age equation for all adults.",
  "  Ya(%) = 99.06 - 0.02*Age.",
  "",
  "CrapoHybrid:",
  "  Pediatric (<18 y): measured individual Ya when available; otherwise Crapo equation extrapolation.",
  "  Adult (>=18 y): Crapo age equation for all adults.",
  "  SaO2(%) = 97.66 - 0.0296*Age.",
  "",
  "The compact input stores precomputed OEF for the primary and two hybrid Ya assumptions.",
  "",
  "Both Ya-derived OEF variables use ORIGINAL Yv so the Ya assumption is tested independently of Hct sensitivity.",
  "Overall/female/male trajectory and growth-rate panels are generated by Step19 from the saved curve tables.",
  "df(mu) is BIC-selected using Original OEF and then locked for both sensitivity models."
)

writeLines(notes, file.path(OutDir, "README_YaSensitivity_Hybrid.txt"))

cat("\n============================================================\n")
cat("Completed. Results saved to:\n", OutDir, "\n")
cat("Run Step19 after Step8 for the overall/female/male figures.\n")
cat("============================================================\n")
