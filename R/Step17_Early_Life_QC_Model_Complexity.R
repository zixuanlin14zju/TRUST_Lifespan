# Step17_Early_Life_QC_Model_Complexity.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 17
# Scientific and validation caveats: docs/ANALYSIS_NOTES.md
source("R/lib/oef_utils.R")
source("R/lib/qc_flags.R")
oef_start_step("Step17")

# ============================================================
# Early-life OEF sensitivity analysis
#
# Primary QC: dR2 <= 10 s^-1
# Strict QC : dR2 <= 5 s^-1

# Main goals:
#   1. Compare early-life trajectories under standard vs strict QC
#   2. Test whether apparent neonatal elevation is QC-dependent
#   3. Test whether an early-life downturn is robust to spline complexity
# ============================================================


# ============================================================
# 00) Packages
# ============================================================

required_packages <- c(
  "readxl",
  "gamlss",
  "gamlss.dist",
  "splines",
  "dplyr",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  stop(
    "Please install the following packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(gamlss)
  library(gamlss.dist)
  library(splines)
  library(dplyr)
  library(ggplot2)
})


# ============================================================
# 01) User settings
# ============================================================

DataDir <- oef_input()

HCFile <- oef_input("trust_hc_hct_ya_sensitivity.xlsx")

OutDir <- file.path(
  oef_output(),
  "EarlyLife_QC_ModelComplexity"
)

if (!dir.exists(OutDir)) {
  dir.create(
    OutDir,
    recursive = TRUE
  )
}

cat(
  "\nResults will be saved to:\n",
  OutDir,
  "\n\n"
)

# ------------------------------------------------------------
# Actual effective cubic B-spline degrees of freedom
# ------------------------------------------------------------

CandidateDfMu <- 3:6

# Final selected sigma complexity
DfSigma <- 3

# GAMLSS control
con <- gamlss.control(
  n.cyc = 200,
  trace = FALSE
)


# ============================================================
# 02) Load FINAL HC dataset
# ============================================================

raw <- read_excel(
  HCFile,
  sheet = "Sheet1"
)

if (!("Sex" %in% names(raw)) && "Sex(M0F1)" %in% names(raw)) {
  raw$Sex <- raw[["Sex(M0F1)"]]
}

required_cols <- c(
  "Age",
  "Sex",
  "SiteID",
  "OEF",
  "dR2_gt5",
  "dR2_gt10"
)

missing_cols <- setdiff(
  required_cols,
  names(raw)
)

if (length(missing_cols) > 0) {
  
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}


# ------------------------------------------------------------
# Construct variables
# ------------------------------------------------------------

qc_flags <- oef_validate_dr2_flags(raw)
raw$dR2_gt5 <- qc_flags$dR2_gt5
raw$dR2_gt10 <- qc_flags$dR2_gt10

dat_all <- raw %>%
  transmute(
    Age = as.numeric(Age),
    
    Sex = as.numeric(
      Sex
    ),
    
    SiteID = factor(
      SiteID
    ),
    
    phenotype = as.numeric(
      OEF
    ),
    
    dR2_gt5 = dR2_gt5,
    dR2_gt10 = dR2_gt10
  ) %>%
  filter(
    !is.na(Age),
    !is.na(Sex),
    !is.na(SiteID),
    !is.na(phenotype),
    !is.na(dR2_gt5),
    !is.na(dR2_gt10)
  ) %>%
  droplevels()


# Check sex coding
if (!all(
  unique(dat_all$Sex) %in% c(0, 1)
)) {
  
  stop(
    "Sex is not cleanly coded as 0/1."
  )
}


cat(
  "Available HC N:",
  nrow(dat_all),
  "\n"
)

cat(
  "Age range:",
  paste(
    round(
      range(dat_all$Age),
      3
    ),
    collapse = " - "
  ),
  "\n"
)

cat(
  "Number of sites:",
  nlevels(dat_all$SiteID),
  "\n\n"
)


# ============================================================
# 03) Define QC datasets
# ============================================================

dat_primary <- dat_all %>%
  filter(
    dR2_gt10 == 0L
  ) %>%
  droplevels()

dat_strict <- dat_all %>%
  filter(
    dR2_gt5 == 0L
  ) %>%
  droplevels()


cat(
  "Primary QC (dR2 <= 10): N =",
  nrow(dat_primary),
  "\n"
)

cat(
  "Strict QC  (dR2 <= 5):  N =",
  nrow(dat_strict),
  "\n\n"
)


cat(
  "Primary QC sites:\n"
)

print(
  table(dat_primary$SiteID)
)

cat(
  "\nStrict QC sites:\n"
)

print(
  table(dat_strict$SiteID)
)


# ============================================================
# 04) Fit one GAMLSS model
# ============================================================

fit_one_model <- function(
    dat,
    df_mu
) {
  
  # Critical:
  # remove unused factor levels independently
  # for each QC dataset
  dat <- droplevels(dat)
  
  fit <- tryCatch(
    
    oef_fit_gamlss(
      
      phenotype ~
        bs(
          Age,
          df = df_mu,
          degree = 3
        ) * Sex +
        random(SiteID),
      
      sigma.fo =
        ~ bs(
          Age,
          df = DfSigma,
          degree = 3
        ) + Sex,
      
      nu.fo =
        ~ 1,
      
      tau.fo =
        ~ 1,
      
      family =
        BCTo,
      
      data =
        dat,
      
      control =
        con
    ),
    
    error = function(e) {
      
      message(
        "Model fitting failed at df_mu = ",
        df_mu,
        ": ",
        e$message
      )
      
      NULL
    }
  )
  
  return(
    fit
  )
}


# ============================================================
# 05) Fit full df set for one QC condition
# ============================================================

fit_model_set <- function(
    dat,
    QC_label
) {
  
  dat <- droplevels(dat)
  
  cat(
    "\n============================================\n"
  )
  
  cat(
    QC_label,
    "\n"
  )
  
  cat(
    "============================================\n"
  )
  
  
  models <- vector(
    "list",
    length(CandidateDfMu)
  )
  
  names(models) <- paste0(
    "df",
    CandidateDfMu
  )
  
  
  result_list <- vector(
    "list",
    length(CandidateDfMu)
  )
  
  
  for (i in seq_along(
    CandidateDfMu
  )) {
    
    df_mu <- CandidateDfMu[i]
    
    cat(
      "\nFitting mu df =",
      df_mu,
      "...\n"
    )
    
    mod <- fit_one_model(
      dat = dat,
      df_mu = df_mu
    )
    
    models[i] <- list(mod)
    
    
    if (is.null(mod)) {
      
      result_list[[i]] <- data.frame(
        
        QC =
          QC_label,
        
        DF_Mu =
          df_mu,
        
        DF_Sigma =
          DfSigma,
        
        N =
          nrow(dat),
        
        N_Sites =
          nlevels(dat$SiteID),
        
        Converged =
          FALSE,
        
        GlobalDeviance =
          NA_real_,
        
        AIC =
          NA_real_,
        
        BIC =
          NA_real_
      )
      
    } else {
      
      result_list[[i]] <- data.frame(
        
        QC =
          QC_label,
        
        DF_Mu =
          df_mu,
        
        DF_Sigma =
          DfSigma,
        
        N =
          nrow(dat),
        
        N_Sites =
          nlevels(dat$SiteID),
        
        Converged =
          isTRUE(
            mod$converged
          ),
        
        GlobalDeviance =
          mod$G.deviance,
        
        AIC =
          AIC(mod),
        
        BIC =
          mod$sbc
      )
    }
  }
  
  
  results_df <- bind_rows(
    result_list
  )
  
  
  valid_bic <- results_df$BIC[
    results_df$Converged & is.finite(results_df$BIC)
  ]
  
  
  if (length(valid_bic) > 0) {
    
    results_df$Delta_BIC <-
      results_df$BIC -
      min(valid_bic)
    
  } else {
    
    results_df$Delta_BIC <-
      NA_real_
  }
  
  
  cat(
    "\nModel comparison:\n"
  )
  
  print(
    results_df
  )
  
  
  list(
    models = models,
    results = results_df,
    data = dat
  )
}


# ============================================================
# 06) Run both QC analyses
# ============================================================

fit_primary <- fit_model_set(
  dat = dat_primary,
  QC_label =
    "Primary QC (dR2 <= 10)"
)

fit_strict <- fit_model_set(
  dat = dat_strict,
  QC_label =
    "Strict QC (dR2 <= 5)"
)


model_compare <- bind_rows(
  fit_primary$results,
  fit_strict$results
)


write.csv(
  model_compare,
  file.path(
    OutDir,
    "01_model_comparison_primary_vs_strict_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 07) Best model under each QC
# ============================================================

best_models <- model_compare %>%
  filter(
    Converged,
    is.finite(BIC)
  ) %>%
  group_by(
    QC
  ) %>%
  slice_min(
    order_by = BIC,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()


cat(
  "\n============================================\n"
)

cat(
  "BEST MODEL BY BIC\n"
)

cat(
  "============================================\n"
)

print(
  best_models
)


write.csv(
  best_models,
  file.path(
    OutDir,
    "02_best_model_by_QC.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 08) Early-life descriptive statistics
# Continuous delta-R2 mean/SD cannot be calculated from threshold flags.
# ============================================================

make_early_summary <- function(
    dat,
    QC_label
) {
  
  dat <- droplevels(dat)
  
  
  out <- dat %>%
    
    filter(
      Age <= 10
    ) %>%
    
    mutate(
      
      AgeGroup = cut(
        
        Age,
        
        breaks = c(
          -Inf,
          1,
          2,
          5,
          10.000001
        ),
        
        right = FALSE,
        
        labels = c(
          "<1",
          "1-<2",
          "2-<5",
          "5-10"
        )
      )
    ) %>%
    
    group_by(
      AgeGroup
    ) %>%
    
    summarise(
      
      QC =
        QC_label,
      
      N =
        n(),
      
      Age_mean =
        mean(
          Age,
          na.rm = TRUE
        ),
      
      Age_SD =
        sd(
          Age,
          na.rm = TRUE
        ),
      
      OEF_mean =
        mean(
          phenotype,
          na.rm = TRUE
        ),
      
      OEF_SD =
        sd(
          phenotype,
          na.rm = TRUE
        ),
      
      OEF_median =
        median(
          phenotype,
          na.rm = TRUE
        ),
      
      Site1_N =
        sum(
          as.character(
            SiteID
          ) == "1"
        ),
      
      Site1_percent =
        100 *
        mean(
          as.character(
            SiteID
          ) == "1"
        ),
      
      .groups =
        "drop"
    )
  
  return(
    out
  )
}


summary_primary <- make_early_summary(
  dat_primary,
  "Primary QC (dR2 <= 10)"
)

summary_strict <- make_early_summary(
  dat_strict,
  "Strict QC (dR2 <= 5)"
)


early_summary <- bind_rows(
  summary_primary,
  summary_strict
)


cat(
  "\n============================================\n"
)

cat(
  "EARLY-LIFE RAW SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  early_summary
)


write.csv(
  early_summary,
  file.path(
    OutDir,
    "03_early_life_raw_summary_primary_vs_strict.csv"
  ),
  row.names = FALSE
)


# ============================================================
# ============================================================

predict_median_safe <- function(
    model,
    age_grid,
    sex_value,
    ref_data
) {
  
  if (is.null(model)) {
    
    return(
      rep(
        NA_real_,
        length(age_grid)
      )
    )
  }
  
  
  ref_data <- droplevels(
    ref_data
  )
  
  
  site_levels_local <-
    levels(
      ref_data$SiteID
    )
  
  
  if (length(site_levels_local) == 0) {
    
    stop(
      "No valid SiteID levels."
    )
  }
  
  
  site_ref <-
    site_levels_local[1]
  
  
  nd <- data.frame(
    
    Age =
      age_grid,
    
    Sex =
      rep(
        sex_value,
        length(age_grid)
      ),
    
    SiteID =
      factor(
        rep(
          site_ref,
          length(age_grid)
        ),
        levels =
          site_levels_local
      )
  )
  
  
  pred <- tryCatch(
    
    oef_predict_all(random = "zero", 
      model,
      newdata = nd
    ),
    
    error = function(e) {
      
      message(
        "Prediction failed: ",
        e$message
      )
      
      NULL
    }
  )
  
  
  if (is.null(pred)) {
    
    return(
      rep(
        NA_real_,
        length(age_grid)
      )
    )
  }
  
  
  med <- tryCatch(
    
    qBCTo(
      
      0.5,
      
      mu =
        pred$mu,
      
      sigma =
        pred$sigma,
      
      nu =
        pred$nu,
      
      tau =
        pred$tau
    ),
    
    error = function(e) {
      
      message(
        "BCTo median calculation failed: ",
        e$message
      )
      
      rep(
        NA_real_,
        length(age_grid)
      )
    }
  )
  
  
  return(
    as.numeric(
      med
    )
  )
}


# ============================================================
# 10) Generate trajectory set
# ============================================================

make_trajectory_set <- function(
    model_list,
    ref_data,
    QC_label
) {
  
  ref_data <- droplevels(
    ref_data
  )
  
  
  age_grid <- seq(
    0,
    10,
    by = 0.02
  )
  
  
  out <- vector(
    "list",
    length(CandidateDfMu)
  )
  
  
  for (i in seq_along(
    CandidateDfMu
  )) {
    
    df_mu <-
      CandidateDfMu[i]
    
    mod <-
      model_list[[i]]
    
    
    cat(
      "Predicting ",
      QC_label,
      ", df_mu = ",
      df_mu,
      "...\n",
      sep = ""
    )
    
    
    pred_male <- predict_median_safe(
      
      model =
        mod,
      
      age_grid =
        age_grid,
      
      sex_value =
        0,
      
      ref_data =
        ref_data
    )
    
    
    pred_female <- predict_median_safe(
      
      model =
        mod,
      
      age_grid =
        age_grid,
      
      sex_value =
        1,
      
      ref_data =
        ref_data
    )
    
    
    pred_avg <-
      (
        pred_male +
          pred_female
      ) / 2
    
    
    out[[i]] <- data.frame(
      
      QC =
        QC_label,
      
      DF_Mu =
        df_mu,
      
      Age =
        age_grid,
      
      Male =
        pred_male,
      
      Female =
        pred_female,
      
      OEF =
        pred_avg
    )
  }
  
  
  bind_rows(
    out
  )
}


# ============================================================
# 11) Generate trajectories
# ============================================================

cat(
  "\nGenerating primary-QC trajectories...\n"
)

traj_primary <- make_trajectory_set(
  
  model_list =
    fit_primary$models,
  
  ref_data =
    fit_primary$data,
  
  QC_label =
    "Primary QC (dR2 <= 10)"
)


cat(
  "\nGenerating strict-QC trajectories...\n"
)

traj_strict <- make_trajectory_set(
  
  model_list =
    fit_strict$models,
  
  ref_data =
    fit_strict$data,
  
  QC_label =
    "Strict QC (dR2 <= 5)"
)


trajectory_all <- bind_rows(
  traj_primary,
  traj_strict
)


write.csv(
  trajectory_all,
  file.path(
    OutDir,
    "04_early_life_trajectories_primary_vs_strict.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12) Check prediction completeness
# ============================================================

prediction_check <- trajectory_all %>%
  
  group_by(
    QC,
    DF_Mu
  ) %>%
  
  summarise(
    
    N_predictions =
      n(),
    
    N_valid =
      sum(
        is.finite(OEF)
      ),
    
    N_missing =
      sum(
        !is.finite(OEF)
      ),
    
    .groups =
      "drop"
  )


cat(
  "\n============================================\n"
)

cat(
  "PREDICTION CHECK\n"
)

cat(
  "============================================\n"
)

print(
  prediction_check
)


write.csv(
  prediction_check,
  file.path(
    OutDir,
    "05_prediction_check.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 13) Raw early-life data for plotting
# ============================================================

raw_primary_plot <- dat_primary %>%
  
  filter(
    Age <= 10
  ) %>%
  
  mutate(
    QC =
      "Primary QC (dR2 <= 10)"
  )


raw_strict_plot <- dat_strict %>%
  
  filter(
    Age <= 10
  ) %>%
  
  mutate(
    QC =
      "Strict QC (dR2 <= 5)"
  )


raw_plot <- bind_rows(
  raw_primary_plot,
  raw_strict_plot
)


trajectory_all <- trajectory_all %>%
  
  mutate(
    
    DF_Label = factor(
      
      paste0(
        "df = ",
        DF_Mu
      ),
      
      levels = paste0(
        "df = ",
        CandidateDfMu
      )
    ),
    
    QC = factor(
      
      QC,
      
      levels = c(
        "Primary QC (dR2 <= 10)",
        "Strict QC (dR2 <= 5)"
      )
    )
  )


raw_plot$QC <- factor(
  
  raw_plot$QC,
  
  levels = c(
    "Primary QC (dR2 <= 10)",
    "Strict QC (dR2 <= 5)"
  )
)


# ============================================================
# 14) Plot A:
# df = 3-6 under both QC thresholds
# ============================================================

p_all <- ggplot() +
  
  geom_point(
    
    data =
      raw_plot,
    
    aes(
      x = Age,
      y = phenotype
    ),
    
    alpha =
      0.28,
    
    size =
      1.6
  ) +
  
  geom_line(
    
    data =
      trajectory_all %>%
      filter(
        is.finite(OEF)
      ),
    
    aes(
      x = Age,
      y = OEF,
      linetype = DF_Label
    ),
    
    linewidth =
      1.0
  ) +
  
  facet_wrap(
    ~ QC,
    ncol = 1
  ) +
  
  labs(
    
    x =
      "Age (years)",
    
    y =
      "OEF",
    
    linetype =
      expression(df[mu]),
    
    title =
      "Early-life OEF trajectory by model complexity and QC threshold"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      ),
    
    strip.text =
      element_text(
        face = "bold"
      )
  )


print(
  p_all
)


ggsave(
  
  filename = file.path(
    OutDir,
    "06_early_life_all_df_primary_vs_strict_QC.png"
  ),
  
  plot =
    p_all,
  
  width =
    8,
  
  height =
    9,
  
  dpi =
    300
)


# ============================================================
# 15) Plot B:
# Direct comparison of primary df_mu = 3
# ============================================================

traj_df3 <- trajectory_all %>%
  
  filter(
    DF_Mu == 3,
    is.finite(OEF)
  )


p_df3 <- ggplot() +
  
  geom_point(
    
    data =
      dat_primary %>%
      filter(
        Age <= 10
      ),
    
    aes(
      x = Age,
      y = phenotype
    ),
    
    alpha =
      0.20,
    
    size =
      1.5
  ) +
  
  geom_line(
    
    data =
      traj_df3,
    
    aes(
      x = Age,
      y = OEF,
      linetype = QC
    ),
    
    linewidth =
      1.2
  ) +
  
  labs(
    
    x =
      "Age (years)",
    
    y =
      "OEF",
    
    linetype =
      "QC",
    
    title =
      "Early-life OEF trajectory under primary and strict QC"
  ) +
  
  theme_classic(
    base_size = 14
  ) +
  
  theme(
    
    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      )
  )


print(
  p_df3
)


ggsave(
  
  filename = file.path(
    OutDir,
    "07_df3_primary_vs_strict_QC_overlay.png"
  ),
  
  plot =
    p_df3,
  
  width =
    8,
  
  height =
    6,
  
  dpi =
    300
)


# ============================================================
# 16) Quantify trajectory differences
# ============================================================

traj_difference <- traj_primary %>%
  
  select(
    DF_Mu,
    Age,
    OEF_Primary = OEF
  ) %>%
  
  inner_join(
    
    traj_strict %>%
      
      select(
        DF_Mu,
        Age,
        OEF_Strict = OEF
      ),
    
    by = c(
      "DF_Mu",
      "Age"
    )
  ) %>%
  
  mutate(
    
    Difference_Strict_minus_Primary =
      OEF_Strict -
      OEF_Primary,
    
    Abs_Difference =
      abs(
        Difference_Strict_minus_Primary
      )
  )


write.csv(
  
  traj_difference,
  
  file.path(
    OutDir,
    "08_trajectory_difference_strict_minus_primary.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 17) Summarize trajectory differences
# ============================================================

trajectory_difference_summary <-
  
  traj_difference %>%
  
  group_by(
    DF_Mu
  ) %>%
  
  summarise(
    
    Mean_difference =
      mean(
        Difference_Strict_minus_Primary,
        na.rm = TRUE
      ),
    
    Mean_absolute_difference =
      mean(
        Abs_Difference,
        na.rm = TRUE
      ),
    
    Maximum_absolute_difference =
      max(
        Abs_Difference,
        na.rm = TRUE
      ),
    
    Difference_at_age0 =
      Difference_Strict_minus_Primary[
        which.min(
          abs(
            Age - 0
          )
        )
      ],
    
    Difference_at_age1 =
      Difference_Strict_minus_Primary[
        which.min(
          abs(
            Age - 1
          )
        )
      ],
    
    Difference_at_age5 =
      Difference_Strict_minus_Primary[
        which.min(
          abs(
            Age - 5
          )
        )
      ],
    
    Difference_at_age10 =
      Difference_Strict_minus_Primary[
        which.min(
          abs(
            Age - 10
          )
        )
      ],
    
    .groups =
      "drop"
  )


cat(
  "\n============================================\n"
)

cat(
  "TRAJECTORY DIFFERENCE SUMMARY\n"
)

cat(
  "============================================\n"
)

print(
  trajectory_difference_summary
)


write.csv(
  
  trajectory_difference_summary,
  
  file.path(
    OutDir,
    "09_trajectory_difference_summary.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 18) Simple sample-count comparison
# ============================================================

sample_counts <- data.frame(
  
  QC = c(
    "Primary QC (dR2 <= 10)",
    "Strict QC (dR2 <= 5)"
  ),
  
  Total_N = c(
    nrow(dat_primary),
    nrow(dat_strict)
  ),
  
  N_Age_le10 = c(
    
    sum(
      dat_primary$Age <= 10
    ),
    
    sum(
      dat_strict$Age <= 10
    )
  ),
  
  N_Age_lt1 = c(
    
    sum(
      dat_primary$Age < 1
    ),
    
    sum(
      dat_strict$Age < 1
    )
  ),
  
  N_Age_1_to_lt2 = c(
    
    sum(
      dat_primary$Age >= 1 &
        dat_primary$Age < 2
    ),
    
    sum(
      dat_strict$Age >= 1 &
        dat_strict$Age < 2
    )
  ),
  
  N_Age_2_to_lt5 = c(
    
    sum(
      dat_primary$Age >= 2 &
        dat_primary$Age < 5
    ),
    
    sum(
      dat_strict$Age >= 2 &
        dat_strict$Age < 5
    )
  ),
  
  N_Age_5_to_10 = c(
    
    sum(
      dat_primary$Age >= 5 &
        dat_primary$Age <= 10
    ),
    
    sum(
      dat_strict$Age >= 5 &
        dat_strict$Age <= 10
    )
  )
)


cat(
  "\n============================================\n"
)

cat(
  "SAMPLE COUNTS\n"
)

cat(
  "============================================\n"
)

print(
  sample_counts
)


write.csv(
  
  sample_counts,
  
  file.path(
    OutDir,
    "10_sample_counts_by_QC.csv"
  ),
  
  row.names = FALSE
)


# ============================================================
# 19) Save fitted models
# ============================================================

for (i in seq_along(
  CandidateDfMu
)) {
  
  df_mu <-
    CandidateDfMu[i]
  
  
  if (!is.null(
    fit_primary$models[[i]]
  )) {
    
    saveRDS(
      
      fit_primary$models[[i]],
      
      file = file.path(
        OutDir,
        paste0(
          "model_primaryQC_dfmu",
          df_mu,
          ".rds"
        )
      )
    )
  }
  
  
  if (!is.null(
    fit_strict$models[[i]]
  )) {
    
    saveRDS(
      
      fit_strict$models[[i]],
      
      file = file.path(
        OutDir,
        paste0(
          "model_strictQC_dfmu",
          df_mu,
          ".rds"
        )
      )
    )
  }
}


# ============================================================
# 20) Finish
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "ANALYSIS COMPLETED\n"
)

cat(
  "============================================\n"
)

cat(
  "Results saved to:\n",
  OutDir,
  "\n"
)
