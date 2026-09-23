# Step10_Hct_Sensitivity_Disease.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 10
source("R/lib/oef_utils.R")
oef_start_step("Step10")

# ============================================================
# TRUST lifespan OEF Hct sensitivity:
# patient Z-scores and disease-pattern analysis
#
# This script intentionally mirrors the ORIGINAL Step3 workflow:
#   1) for each Hct scenario, refit a clean BCTo GAMLSS model on HC
#      with df_mu = 3 and df_sigma = 3;
#   2) predict patient distribution parameters;
#   3) convert observed patient OEF -> BCTo CDF -> normal Z-score;
#   4) retain diagnoses with N >= 5;
#   5) calculate descriptive Z statistics;
#   6) for each diagnosis run:
#        - one-sample t-test vs 0
#        - z ~ Age + Sex + (1|SiteID) when >1 site
#          otherwise z ~ Age + Sex
#
# Four requested Hct sensitivity scenarios:
#   RegionMatched, Zierk, MahlknechtMean, CohortEmpirical
#
# Additional sensitivity-comparison tables are produced after the
# Step3-equivalent analysis; these do NOT change the Step3 statistics.
# ============================================================

# ============================================================
# 00) Packages
# ============================================================

required_packages <- c(
  "gamlss",
  "dplyr",
  "lme4",
  "lmerTest",
  "splines",
  "readxl"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(gamlss)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(splines)
  library(readxl)
})

# ============================================================
# 01) User settings
# ============================================================

DataDir <- oef_input()

HCPath <- oef_input("trust_hc_hct_ya_sensitivity.xlsx")

PatientPath <- oef_input("trust_patient_hct_sensitivity.xlsx")

OutDir <- file.path(
  oef_output(),
  "GAMLSS_HctSensitivity_Disease"
)

if (!dir.exists(OutDir)) {
  dir.create(OutDir, recursive = TRUE)
}

DiagnosisCol <- "Diagnosis"
MinSamplesPerDx <- 5
Zthr <- 1.96

# EXACTLY as in original Step3
df_mu <- 3
df_sigma <- 3
con <- gamlss.control(n.cyc = 200, trace = FALSE)

ScenarioMap <- c(
  Original        = "OEF",
  RegionMatched   = "OEF_RegionMatched",
  Zierk           = "OEF_Zierk",
  MahlknechtMean  = "OEF_MahlMean",
  CohortEmpirical = "OEF_CohortEmpirical"
)

SensitivityScenarios <- setdiff(
  names(ScenarioMap),
  "Original"
)

# ============================================================
# 02) Helpers
# ============================================================

normalize_sex <- function(dat) {

  if ("Sex(M0F1)" %in% names(dat)) {
    x <- dat[["Sex(M0F1)"]]
  } else if ("Sex" %in% names(dat)) {
    x <- dat[["Sex"]]
  } else {
    stop("No Sex or Sex(M0F1) column found.")
  }

  if (is.numeric(x) || is.integer(x)) {
    out <- as.numeric(x)

    if (!all(na.omit(unique(out)) %in% c(0, 1))) {
      stop("Numeric Sex is not coded 0/1.")
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
    stop("Could not map all Sex values.")
  }

  out
}

read_current_excel <- function(path, label) {

  if (!file.exists(path)) {
    stop(label, " file not found: ", path)
  }

  dat <- as.data.frame(
    readxl::read_excel(path, sheet = "Sheet1"),
    stringsAsFactors = FALSE
  )

  dat$Sex <- normalize_sex(dat)
  dat$Age <- as.numeric(dat$Age)

  if (!("SiteID" %in% names(dat))) {
    stop(label, ": SiteID column missing.")
  }

  if ("QC" %in% names(dat)) {
    keep <- is.na(dat$QC) |
      toupper(trimws(as.character(dat$QC))) == "PASS"

    cat(
      sprintf(
        "%s QC: retained %d/%d rows\n",
        label,
        sum(keep),
        nrow(dat)
      )
    )

    dat <- dat[keep, , drop = FALSE]
  }

  dat$.RowID <- seq_len(nrow(dat))

  dat
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok])
}

safe_mae <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  mean(abs(x[ok] - y[ok]))
}

# ============================================================
# 03) Load data
# ============================================================

M_HC <- read_current_excel(HCPath, "HC")
patient_master <- read_current_excel(PatientPath, "Patient")

required_hc <- c(
  "Age",
  "Sex",
  "SiteID",
  unname(ScenarioMap)
)

required_pt <- c(
  "Age",
  "Sex",
  "SiteID",
  DiagnosisCol,
  unname(ScenarioMap)
)

missing_hc <- setdiff(required_hc, names(M_HC))
missing_pt <- setdiff(required_pt, names(patient_master))

if (length(missing_hc) > 0) {
  stop(
    "HC missing: ",
    paste(missing_hc, collapse = ", ")
  )
}

if (length(missing_pt) > 0) {
  stop(
    "Patient missing: ",
    paste(missing_pt, collapse = ", ")
  )
}

cat("\n=== Patient data overview ===\n")
cat(
  sprintf(
    "Age range: %.1f - %.1f years\n",
    min(patient_master$Age, na.rm = TRUE),
    max(patient_master$Age, na.rm = TRUE)
  )
)

cat(
  sprintf(
    "Sex counts: F=%d, M=%d\n",
    sum(patient_master$Sex == 1, na.rm = TRUE),
    sum(patient_master$Sex == 0, na.rm = TRUE)
  )
)

cat("Diagnosis counts:\n")
print(table(patient_master[[DiagnosisCol]]))

# ============================================================
# 04) Run original Step3-equivalent analysis for each scenario
# ============================================================

AllPatientZ <- patient_master
AllDesc <- list()
AllResults <- list()
ModelList <- list()

for (scenario in names(ScenarioMap)) {

  phenotype_col <- ScenarioMap[[scenario]]

  cat("\n\n")
  cat("============================================================\n")
  cat("SCENARIO:", scenario, "\n")
  cat("Patient phenotype:", phenotype_col, "\n")
  cat("============================================================\n")

  # ----------------------------------------------------------
  # 4.1 Scenario-specific HC dataset
  # ----------------------------------------------------------

  hc <- M_HC
  hc$phenotype <- as.numeric(hc[[phenotype_col]])

  hc <- hc[, c("Age", "Sex", "SiteID", "phenotype"), drop = FALSE]

  hc <- hc[
    is.finite(hc$Age) &
    !is.na(hc$Sex) &
    !is.na(hc$SiteID) &
    is.finite(hc$phenotype),
    ,
    drop = FALSE
  ]

  hc <- na.omit(hc)

  # This is equivalent to levels(factor(M_HC$SiteID)) in original Step3.
  hc$SiteID <- factor(hc$SiteID)
  train_site_levels <- levels(hc$SiteID)

  # ----------------------------------------------------------
  # 4.2 Patient data / SiteID handling EXACTLY following Step3
  # ----------------------------------------------------------

  patient_data <- patient_master
  patient_data$phenotype <- as.numeric(
    patient_data[[phenotype_col]]
  )

  patient_data$SiteID <- factor(
    patient_data$SiteID,
    levels = train_site_levels
  )

  newdata_pat <- data.frame(
    Age = patient_data$Age,
    Sex = patient_data$Sex
  )

  site_dummy <- train_site_levels[
    which(!is.na(train_site_levels))[1]
  ]

  newdata_pat$SiteID <- factor(
    rep(site_dummy, nrow(newdata_pat)),
    levels = train_site_levels
  )

  # ----------------------------------------------------------
  # 4.3 Refit clean GAMLSS on HC EXACTLY as original Step3
  # ----------------------------------------------------------

  cat("\n=== Refitting clean GAMLSS model on HC ===\n")

  clean_model <- oef_fit_gamlss(
    phenotype ~
      bs(Age, df = df_mu) * Sex +
      random(as.factor(SiteID)),
    sigma.fo = ~ bs(Age, df = df_sigma) + Sex,
    nu.fo = ~ 1,
    tau.fo = ~ 1,
    family = BCTo,
    data = hc,
    control = con
  )

  ModelList[[scenario]] <- clean_model

  saveRDS(
    clean_model,
    file.path(
      OutDir,
      paste0(
        "clean_model_",
        scenario,
        ".rds"
      )
    )
  )

  # ----------------------------------------------------------
  # 4.4 Patient Z-score EXACTLY as original Step3
  # ----------------------------------------------------------

  cat(
    "\n=== Predicting distribution parameters ===\n"
  )

  pred_all <- oef_predict_all(random = "zero", 
    clean_model,
    newdata = newdata_pat
  )

  mu_patient <- pred_all$mu
  sigma_patient <- pred_all$sigma
  nu_patient <- pred_all$nu
  tau_patient <- pred_all$tau

  y <- patient_data$phenotype

  cumulative_prob <- pBCTo(
    y,
    mu = mu_patient,
    sigma = sigma_patient,
    nu = nu_patient,
    tau = tau_patient
  )

  cumulative_prob <- pmin(
    pmax(cumulative_prob, 1e-10),
    1 - 1e-10
  )

  patient_data$z_score <- qnorm(
    cumulative_prob
  )

  AllPatientZ[[paste0("z_", scenario)]] <- patient_data$z_score

  cat("\n=== Patient Z-score summary ===\n")
  cat(
    sprintf(
      "Mean: %.3f | SD: %.3f | Range: [%.3f, %.3f]\n",
      mean(patient_data$z_score, na.rm = TRUE),
      sd(patient_data$z_score, na.rm = TRUE),
      min(patient_data$z_score, na.rm = TRUE),
      max(patient_data$z_score, na.rm = TRUE)
    )
  )

  # Save scenario-specific patient Z file
  write.csv(
    patient_data,
    file.path(
      OutDir,
      paste0(
        "disease_zscores_",
        scenario,
        "_reference_site.csv"
      )
    ),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # 4.5 Diagnosis filtering EXACTLY as original Step3
  # ----------------------------------------------------------

  dx_counts <- table(
    patient_data[[DiagnosisCol]]
  )

  valid_dx <- names(
    dx_counts[
      dx_counts >= MinSamplesPerDx
    ]
  )

  cat("\n=== Diagnosis filtering ===\n")
  cat(
    sprintf(
      "Total diagnosis categories: %d\n",
      length(dx_counts)
    )
  )
  cat(
    sprintf(
      "Keeping dx with N >= %d: %d categories\n",
      MinSamplesPerDx,
      length(valid_dx)
    )
  )

  patient_f <- patient_data %>%
    filter(
      .data[[DiagnosisCol]] %in%
        valid_dx
    ) %>%
    mutate(
      Diagnosis = factor(
        .data[[DiagnosisCol]],
        levels = valid_dx
      )
    )

  if (length(valid_dx) < 2) {
    stop(
      "Need at least 2 diagnosis groups after filtering."
    )
  }

  # ----------------------------------------------------------
  # 4.6 Descriptive statistics EXACTLY as original Step3
  # ----------------------------------------------------------

  desc_stats <- patient_f %>%
    group_by(Diagnosis) %>%
    summarise(
      N = n(),
      Mean_Z = mean(
        z_score,
        na.rm = TRUE
      ),
      SD_Z = sd(
        z_score,
        na.rm = TRUE
      ),
      SE_Z = SD_Z / sqrt(N),
      Median_Z = median(
        z_score,
        na.rm = TRUE
      ),
      Min_Z = min(
        z_score,
        na.rm = TRUE
      ),
      Max_Z = max(
        z_score,
        na.rm = TRUE
      ),
      CI95_L = Mean_Z -
        Zthr * SE_Z,
      CI95_U = Mean_Z +
        Zthr * SE_Z,
      Abnormal_percent =
        mean(
          abs(z_score) > Zthr,
          na.rm = TRUE
        ) * 100,
      .groups = "drop"
    ) %>%
    arrange(desc(Mean_Z))

  desc_stats$Scenario <- scenario

  AllDesc[[scenario]] <- desc_stats

  write.csv(
    desc_stats,
    file.path(
      OutDir,
      paste0(
        "disease_descriptive_",
        scenario,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  # ----------------------------------------------------------
  # 4.7 Statistics EXACTLY as original Step3
  # ----------------------------------------------------------

  results_list <- list()

  for (
    dx in levels(
      patient_f$Diagnosis
    )
  ) {

    diag_data <- patient_f %>%
      filter(Diagnosis == dx)

    n_sites <- if (
      "SiteID" %in% names(diag_data)
    ) {
      length(
        unique(
          diag_data$SiteID[
            !is.na(diag_data$SiteID)
          ]
        )
      )
    } else {
      0
    }

    cat(
      "\n----------------------------------------\n"
    )

    cat(
      sprintf(
        "Diagnosis: %s | N=%d | Sites=%d\n",
        dx,
        nrow(diag_data),
        n_sites
      )
    )

    # Simple one-sample t-test
    t_res <- t.test(
      diag_data$z_score,
      mu = 0
    )

    # Same fixed-effects construction as original Step3.
    base_terms <- c("1")

    if ("Age" %in% names(diag_data)) {
      base_terms <- c(
        base_terms,
        "Age"
      )
    }

    if ("Sex" %in% names(diag_data)) {
      base_terms <- c(
        base_terms,
        "Sex"
      )
    }

    fixed_part <- paste(
      base_terms,
      collapse = " + "
    )

    use_mixed <-
      ("SiteID" %in% names(diag_data)) &&
      (n_sites > 1)

    model_type <- NA_character_
    model_intercept <- NA_real_
    model_se <- NA_real_
    model_t <- NA_real_
    model_p <- NA_real_
    age_beta <- NA_real_
    age_p <- NA_real_
    sex_beta <- NA_real_
    sex_p <- NA_real_

    if (use_mixed) {

      model_type <-
        "Mixed model: z ~ Age + Sex + (1|SiteID)"

      fml <- as.formula(
        paste0(
          "z_score ~ ",
          fixed_part,
          " + (1|SiteID)"
        )
      )

      mixed_ok <- TRUE

      tryCatch({

        m <- lmer(
          fml,
          data = diag_data
        )

        ct <- coef(summary(m))

        model_intercept <-
          ct[
            "(Intercept)",
            "Estimate"
          ]

        model_se <-
          ct[
            "(Intercept)",
            "Std. Error"
          ]

        model_t <-
          ct[
            "(Intercept)",
            "t value"
          ]

        model_p <-
          ct[
            "(Intercept)",
            "Pr(>|t|)"
          ]

        if ("Age" %in% rownames(ct)) {
          age_beta <-
            ct["Age", "Estimate"]
          age_p <-
            ct["Age", "Pr(>|t|)"]
        }

        sex_rows <- rownames(ct)[
          grepl("^Sex", rownames(ct))
        ]

        if (length(sex_rows) > 0) {
          sex_beta <-
            ct[
              sex_rows[1],
              "Estimate"
            ]
          sex_p <-
            ct[
              sex_rows[1],
              "Pr(>|t|)"
            ]
        }

      }, error = function(e) {

        cat(
          sprintf(
            "Mixed model failed: %s\n",
            e$message
          )
        )

        mixed_ok <<- FALSE
      })

      if (!mixed_ok) {
        use_mixed <- FALSE
      }
    }

    if (!use_mixed) {

      model_type <-
        "Linear regression: z ~ Age + Sex"

      fml <- as.formula(
        paste0(
          "z_score ~ ",
          fixed_part
        )
      )

      m <- lm(
        fml,
        data = diag_data
      )

      ct <- coef(summary(m))

      model_intercept <-
        ct[
          "(Intercept)",
          "Estimate"
        ]

      model_se <-
        ct[
          "(Intercept)",
          "Std. Error"
        ]

      model_t <-
        ct[
          "(Intercept)",
          "t value"
        ]

      model_p <-
        ct[
          "(Intercept)",
          "Pr(>|t|)"
        ]

      if ("Age" %in% rownames(ct)) {
        age_beta <-
          ct["Age", "Estimate"]
        age_p <-
          ct["Age", "Pr(>|t|)"]
      }

      sex_rows <- rownames(ct)[
        grepl("^Sex", rownames(ct))
      ]

      if (length(sex_rows) > 0) {
        sex_beta <-
          ct[
            sex_rows[1],
            "Estimate"
          ]
        sex_p <-
          ct[
            sex_rows[1],
            "Pr(>|t|)"
          ]
      }
    }

    results_list[[as.character(dx)]] <- data.frame(
      Scenario = scenario,
      Diagnosis = as.character(dx),
      N = nrow(diag_data),
      N_sites = ifelse(
        "SiteID" %in% names(diag_data),
        n_sites,
        NA
      ),
      Mean_Z_simple = mean(
        diag_data$z_score,
        na.rm = TRUE
      ),
      SD_Z = sd(
        diag_data$z_score,
        na.rm = TRUE
      ),
      T_test_t = unname(
        t_res$statistic
      ),
      T_test_p = t_res$p.value,
      Model_type = model_type,
      Model_intercept = model_intercept,
      Model_SE = model_se,
      Model_t = model_t,
      Model_p = model_p,
      Age_effect = age_beta,
      Age_p = age_p,
      Sex_effect = sex_beta,
      Sex_p = sex_p,
      Abnormal_percent =
        mean(
          abs(diag_data$z_score) >
            Zthr,
          na.rm = TRUE
        ) * 100
    )
  }

  results_df <- do.call(
    rbind,
    results_list
  ) %>%
    as.data.frame()

  rownames(results_df) <- NULL

  AllResults[[scenario]] <- results_df

  write.csv(
    results_df,
    file.path(
      OutDir,
      paste0(
        "disease_effect_results_",
        scenario,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

# ============================================================
# 05) Save combined Step3-equivalent outputs
# ============================================================

CombinedDesc <- bind_rows(
  AllDesc
)

CombinedResults <- bind_rows(
  AllResults
)

write.csv(
  AllPatientZ,
  file.path(
    OutDir,
    "patient_zscores_all_scenarios.csv"
  ),
  row.names = FALSE
)

write.csv(
  CombinedDesc,
  file.path(
    OutDir,
    "disease_descriptive_all_scenarios.csv"
  ),
  row.names = FALSE
)

write.csv(
  CombinedResults,
  file.path(
    OutDir,
    "disease_effect_results_all_scenarios.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 06) EXTRA sensitivity comparisons vs Original
#     (not part of the original Step3 statistics)
# ============================================================

comparison_rows <- list()

for (scenario in SensitivityScenarios) {

  z_original <-
    AllPatientZ$z_Original

  z_sens <-
    AllPatientZ[[paste0("z_", scenario)]]

  comparison_rows[[scenario]] <-
    data.frame(
      Scenario = scenario,
      N_pair = sum(
        is.finite(z_original) &
          is.finite(z_sens)
      ),
      r_vs_Original =
        safe_cor(
          z_sens,
          z_original
        ),
      MAE_vs_Original =
        safe_mae(
          z_sens,
          z_original
        ),
      MeanDelta_z =
        mean(
          z_sens - z_original,
          na.rm = TRUE
        )
    )
}

PatientSensitivity <-
  bind_rows(comparison_rows)

write.csv(
  PatientSensitivity,
  file.path(
    OutDir,
    "patient_zscore_sensitivity_vs_original.csv"
  ),
  row.names = FALSE
)

# Diagnosis-specific scenario - Original Z shift
delta_rows <- list()

for (scenario in SensitivityScenarios) {

  z_sens_name <-
    paste0("z_", scenario)

  tmp <- AllPatientZ %>%
    mutate(
      .Diagnosis =
        as.character(
          .data[[DiagnosisCol]]
        ),
      .Delta =
        .data[[z_sens_name]] -
        z_Original
    ) %>%
    filter(
      !is.na(.Diagnosis),
      nzchar(.Diagnosis),
      is.finite(.Delta)
    ) %>%
    group_by(.Diagnosis) %>%
    summarise(
      N = n(),
      MeanDelta_z =
        mean(.Delta),
      SDDelta_z =
        sd(.Delta),
      MAEDelta_z =
        mean(abs(.Delta)),
      .groups = "drop"
    ) %>%
    mutate(
      Scenario = scenario
    )

  delta_rows[[scenario]] <- tmp
}

DiagnosisDelta <-
  bind_rows(delta_rows) %>%
  rename(
    Diagnosis = .Diagnosis
  ) %>%
  select(
    Scenario,
    Diagnosis,
    N,
    MeanDelta_z,
    SDDelta_z,
    MAEDelta_z
  )

write.csv(
  DiagnosisDelta,
  file.path(
    OutDir,
    "diagnosis_zscore_change_vs_original.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 07) Console summary
# ============================================================

cat("\n\n============================================================\n")
cat("STEP3 HCT SENSITIVITY ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat("\nOutput directory:\n")
cat(OutDir, "\n")

cat("\nPatient z-score sensitivity vs Original:\n")
print(PatientSensitivity)

cat("\nCombined disease-effect table:\n")
print(
  CombinedResults[
    ,
    c(
      "Scenario",
      "Diagnosis",
      "N",
      "Model_type",
      "Model_intercept",
      "Model_SE",
      "Model_p",
      "T_test_p",
      "Abnormal_percent"
    )
  ]
)

cat("\nDone.\n")
