
suppressPackageStartupMessages({
  library(dplyr)
  library(gamlss)
  library(gamlss.dist)
  library(pracma)
  library(splines)
  library(readr)
})

MainDataPath <- "data/trust_lifespan_hc.csv"

OutDir <- "outputs/main_model_and_sensitivity_compare"
if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)

age_max <- 90
sample_age_step <- 1
main_n_cyc <- 200

main_select_df <- TRUE
main_df_candidates <- 2:5
main_mu_df <- 2
main_sigma_df <- 1

SensitivityFiles <- list(
  StricterQC = list(
    curve = "outputs/sensitivity_stricter_qc/StricterQC_yearly_median_curve.csv",
    growth = "outputs/sensitivity_stricter_qc/StricterQC_yearly_growth_curve.csv"
  ),
  Balanced = list(
    curve = "outputs/sensitivity_balanced_resampling/Balanced_yearly_median_curve.csv",
    growth = "outputs/sensitivity_balanced_resampling/Balanced_yearly_growth_curve.csv"
  ),
  SplitHalf = list(
    curve = "outputs/sensitivity_split_half/SplitHalf_yearly_median_curve_first_vs_second.csv",
    growth = "outputs/sensitivity_split_half/SplitHalf_yearly_growth_first_vs_second.csv"
  ),
  Bootstrap = list(
    curve = "outputs/sensitivity_bootstrap/Bootstrap_yearly_median_curve.csv",
    growth = "outputs/sensitivity_bootstrap/Bootstrap_yearly_growth_curve.csv"
  ),
  LOSO = list(
    curve = "outputs/sensitivity_loso/LOSO_yearly_median_curve.csv",
    growth = "outputs/sensitivity_loso/LOSO_yearly_growth_curve.csv"
  )
)

msg <- function(...) cat(paste0(..., "\n"))

sample_yearly <- function(age, y, by = 1) {
  age_year <- seq(ceiling(min(age, na.rm = TRUE)),
                  floor(max(age, na.rm = TRUE)),
                  by = by)
  y_year <- approx(x = age, y = y, xout = age_year, rule = 2)$y
  data.frame(Age = age_year, Value = y_year)
}

curve_to_growth <- function(y, age) {
  pracma::gradient(y, age)
}

safe_cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) {
    return(list(r = NA_real_, p = NA_real_, n = length(x)))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "pearson"))
  list(r = unname(ct$estimate), p = ct$p.value, n = length(x))
}

safe_rmse <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  sqrt(mean((x[ok] - y[ok])^2))
}

safe_mad <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  mean(abs(x[ok] - y[ok]))
}

safe_maxad <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) return(NA_real_)
  max(abs(x[ok] - y[ok]))
}

unpack_pred <- function(x, n_target = NULL) {
  if (is.null(x)) return(NULL)
  
  if (is.atomic(x) && !is.list(x)) {
    x <- as.numeric(x)
    if (!is.null(n_target) && length(x) == 1) x <- rep(x, n_target)
    return(x)
  }
  
  if (is.list(x)) {
    candidate_names <- c("fit", "fitted.values", "pred", "prediction", "y", "data")
    for (nm in candidate_names) {
      if (!is.null(x[[nm]])) {
        val <- unpack_pred(x[[nm]], n_target = n_target)
        if (!is.null(val)) return(val)
      }
    }
    
    if (length(x) == 1) {
      val <- unpack_pred(x[[1]], n_target = n_target)
      if (!is.null(val)) return(val)
    }
    
    val <- tryCatch(as.numeric(unlist(x)), error = function(e) NULL)
    if (!is.null(val)) {
      if (!is.null(n_target) && length(val) == 1) val <- rep(val, n_target)
      return(val)
    }
    return(NULL)
  }
  
  NULL
}

safe_predict_param <- function(mod, what, newdata, ref_data) {
  out <- tryCatch(
    predict(mod, what = what, newdata = newdata, data = ref_data, type = "response"),
    error = function(e) NULL
  )
  out <- unpack_pred(out, n_target = nrow(newdata))
  
  if (is.null(out) || length(out) == 0) {
    out0 <- tryCatch(
      predict(mod, what = what, data = ref_data, type = "response"),
      error = function(e) NULL
    )
    out0 <- unpack_pred(out0, n_target = nrow(newdata))
    
    if (!is.null(out0) && length(out0) >= 1) {
      return(rep(out0[1], nrow(newdata)))
    } else {
      stop(paste("Failed to predict parameter:", what))
    }
  }
  
  if (length(out) == 1) out <- rep(out, nrow(newdata))
  
  if (length(out) != nrow(newdata)) {
    if (length(out) > nrow(newdata)) {
      out <- out[seq_len(nrow(newdata))]
    } else if (length(out) >= 1) {
      out <- rep(out[1], nrow(newdata))
    } else {
      stop(paste("Predicted parameter has invalid length:", what))
    }
  }
  
  as.numeric(out)
}

predict_median_curve_bcto <- function(mod, age_grid, ref_data, sex_levels = c("0", "1")) {
  ref_data <- ref_data %>%
    dplyr::filter(!is.na(Age), !is.na(Sex), !is.na(SiteID), !is.na(phenotype))
  
  ref_data$SiteID <- factor(ref_data$SiteID)
  ref_data$Sex <- factor(as.character(ref_data$Sex), levels = sex_levels)
  
  site_levels <- levels(ref_data$SiteID)
  if (length(site_levels) == 0) stop("No valid SiteID levels found in ref_data.")
  site_ref <- site_levels[1]
  
  newdata_m <- data.frame(
    Age = age_grid,
    Sex = factor(rep(sex_levels[1], length(age_grid)), levels = sex_levels),
    SiteID = factor(rep(site_ref, length(age_grid)), levels = site_levels)
  )
  
  newdata_f <- data.frame(
    Age = age_grid,
    Sex = factor(rep(sex_levels[2], length(age_grid)), levels = sex_levels),
    SiteID = factor(rep(site_ref, length(age_grid)), levels = site_levels)
  )
  
  mu_m    <- safe_predict_param(mod, "mu",    newdata_m, ref_data)
  sigma_m <- safe_predict_param(mod, "sigma", newdata_m, ref_data)
  nu_m    <- safe_predict_param(mod, "nu",    newdata_m, ref_data)
  tau_m   <- safe_predict_param(mod, "tau",   newdata_m, ref_data)
  
  mu_f    <- safe_predict_param(mod, "mu",    newdata_f, ref_data)
  sigma_f <- safe_predict_param(mod, "sigma", newdata_f, ref_data)
  nu_f    <- safe_predict_param(mod, "nu",    newdata_f, ref_data)
  tau_f   <- safe_predict_param(mod, "tau",   newdata_f, ref_data)
  
  med_m <- qBCTo(0.5, mu = mu_m, sigma = sigma_m, nu = nu_m, tau = tau_m)
  med_f <- qBCTo(0.5, mu = mu_f, sigma = sigma_f, nu = nu_f, tau = tau_f)
  
  list(
    male = med_m,
    female = med_f,
    avg = (med_m + med_f) / 2
  )
}

fit_bcto_model <- function(dat, mu_df = 3, sigma_df = 3, n_cyc = 200, trace = FALSE) {
  dat <- dat %>%
    dplyr::filter(!is.na(Age), !is.na(Sex), !is.na(SiteID), !is.na(phenotype))
  dat$SiteID <- factor(dat$SiteID)
  dat$Sex <- factor(as.character(dat$Sex), levels = c("0", "1"))
  
  gamlss(
    phenotype ~ bs(Age, df = mu_df) * Sex + random(SiteID),
    sigma.fo = ~ bs(Age, df = sigma_df) + Sex,
    nu.fo = ~ 1,
    tau.fo = ~ 1,
    family = BCTo,
    data = dat,
    control = gamlss.control(n.cyc = n_cyc, trace = trace)
  )
}

fit_bcto_select_df <- function(dat, df_candidates = 2:5, sigma_df = 1, n_cyc = 200) {
  dat <- dat %>%
    dplyr::filter(!is.na(Age), !is.na(Sex), !is.na(SiteID), !is.na(phenotype))
  dat$SiteID <- factor(dat$SiteID)
  dat$Sex <- factor(as.character(dat$Sex), levels = c("0", "1"))
  
  models <- lapply(df_candidates, function(df_mu) {
    gamlss(
      phenotype ~ bs(Age, df = df_mu) * Sex + random(SiteID),
      sigma.fo = ~ bs(Age, df = sigma_df) + Sex,
      nu.fo = ~ 1,
      tau.fo = ~ 1,
      family = BCTo,
      data = dat,
      control = gamlss.control(n.cyc = n_cyc, trace = FALSE)
    )
  })
  
  bic_vals <- vapply(models, function(m) m$sbc, numeric(1))
  best_idx <- which.min(bic_vals)
  
  list(
    model = models[[best_idx]],
    best_df_mu = df_candidates[best_idx],
    bic_table = data.frame(df_mu = df_candidates, BIC = bic_vals)
  )
}

read_yearly_curve_file <- function(file, analysis_name, type = c("curve", "growth")) {
  type <- match.arg(type)
  if (!file.exists(file)) stop(paste("File not found:", file))
  
  dat <- read.csv(file, stringsAsFactors = FALSE)
  
  if (all(c("Age", "Value") %in% names(dat))) {
    out <- dat[, c("Age", "Value")]
    names(out)[2] <- "MainLikeValue"
    return(out)
  }
  
  if (all(c("Age", "Value", "Half") %in% names(dat))) {
    dat2 <- dat %>%
      group_by(Age) %>%
      summarise(MainLikeValue = mean(Value, na.rm = TRUE), .groups = "drop")
    return(dat2)
  }
  
  if ("Age" %in% names(dat) && ncol(dat) >= 2) {
    out <- dat[, 1:2]
    names(out) <- c("Age", "MainLikeValue")
    return(out)
  }
  
  stop(paste("Cannot parse file:", file, "for", analysis_name, type))
}

compare_two_yearly_curves <- function(main_df, sens_df, analysis_name, metric_name) {
  merged <- merge(main_df, sens_df, by = "Age", all = FALSE)
  names(merged)[2:3] <- c("Main", "Sensitivity")
  
  ct <- safe_cor_test(merged$Main, merged$Sensitivity)
  
  data.frame(
    Analysis = analysis_name,
    Metric = metric_name,
    N_AgePoints = ct$n,
    r = ct$r,
    p = ct$p,
    rmse = safe_rmse(merged$Main, merged$Sensitivity),
    mean_abs_diff = safe_mad(merged$Main, merged$Sensitivity),
    max_abs_diff = safe_maxad(merged$Main, merged$Sensitivity)
  )
}

msg("=== Load main data ===")
main_dat <- read.csv(MainDataPath, stringsAsFactors = FALSE)

stopifnot(all(c("Age", "Sex", "SiteID", "phenotype") %in% names(main_dat)))

main_dat <- main_dat %>%
  filter(Age <= age_max) %>%
  filter(!is.na(Age), !is.na(Sex), !is.na(SiteID), !is.na(phenotype))

main_dat$SiteID <- factor(main_dat$SiteID)
main_dat$Sex <- as.character(main_dat$Sex)

if (!all(na.omit(unique(main_dat$Sex)) %in% c("0", "1"))) {
  stop("Sex must be coded as 0/1, with 0 = Male and 1 = Female.")
}
main_dat$Sex <- factor(main_dat$Sex, levels = c("0", "1"))

msg(sprintf("N = %d", nrow(main_dat)))
msg(sprintf("Age range = %.2f to %.2f", min(main_dat$Age), max(main_dat$Age)))
msg(sprintf("Sites = %d", nlevels(main_dat$SiteID)))

msg("=== Fit main model ===")
if (main_select_df) {
  main_fit <- fit_bcto_select_df(
    main_dat,
    df_candidates = main_df_candidates,
    sigma_df = main_sigma_df,
    n_cyc = main_n_cyc
  )
  main_model <- main_fit$model
  main_mu_df_used <- main_fit$best_df_mu
  bic_table <- main_fit$bic_table
} else {
  main_model <- fit_bcto_model(
    main_dat,
    mu_df = main_mu_df,
    sigma_df = main_sigma_df,
    n_cyc = main_n_cyc,
    trace = FALSE
  )
  main_mu_df_used <- main_mu_df
  bic_table <- NULL
}

age_grid_dense <- seq(min(main_dat$Age), max(main_dat$Age), by = 0.1)

main_curves <- predict_median_curve_bcto(
  main_model,
  age_grid_dense,
  ref_data = main_dat,
  sex_levels = c("0", "1")
)

main_growth <- curve_to_growth(main_curves$avg, age_grid_dense)
main_growth_f <- curve_to_growth(main_curves$female, age_grid_dense)
main_growth_m <- curve_to_growth(main_curves$male, age_grid_dense)

main_curve_yearly <- sample_yearly(age_grid_dense, main_curves$avg, by = sample_age_step)
main_growth_yearly <- sample_yearly(age_grid_dense, main_growth, by = sample_age_step)

main_curve_yearly_sex <- rbind(
  data.frame(sample_yearly(age_grid_dense, main_curves$female, by = sample_age_step), Sex = "Female"),
  data.frame(sample_yearly(age_grid_dense, main_curves$male,   by = sample_age_step), Sex = "Male")
)

main_growth_yearly_sex <- rbind(
  data.frame(sample_yearly(age_grid_dense, main_growth_f, by = sample_age_step), Sex = "Female"),
  data.frame(sample_yearly(age_grid_dense, main_growth_m, by = sample_age_step), Sex = "Male")
)

peak_idx <- which.max(main_growth)
trough_idx <- which.min(main_growth)

main_summary <- data.frame(
  Analysis = "Main model",
  N_total = nrow(main_dat),
  N_sites = nlevels(main_dat$SiteID),
  Age_min = min(main_dat$Age, na.rm = TRUE),
  Age_max = max(main_dat$Age, na.rm = TRUE),
  Mu_df_used = main_mu_df_used,
  Sigma_df_used = main_sigma_df,
  PeakAge_all = age_grid_dense[peak_idx],
  PeakRate_all = main_growth[peak_idx],
  TroughAge_all = age_grid_dense[trough_idx],
  TroughRate_all = main_growth[trough_idx]
)

write.csv(main_summary,
          file.path(OutDir, "Main_manuscript_summary.csv"),
          row.names = FALSE)

write.csv(main_curve_yearly,
          file.path(OutDir, "Main_yearly_median_curve.csv"),
          row.names = FALSE)

write.csv(main_growth_yearly,
          file.path(OutDir, "Main_yearly_growth_curve.csv"),
          row.names = FALSE)

write.csv(main_curve_yearly_sex,
          file.path(OutDir, "Main_yearly_median_curve_by_sex.csv"),
          row.names = FALSE)

write.csv(main_growth_yearly_sex,
          file.path(OutDir, "Main_yearly_growth_curve_by_sex.csv"),
          row.names = FALSE)

if (!is.null(bic_table)) {
  write.csv(bic_table,
            file.path(OutDir, "Main_model_selection_BIC.csv"),
            row.names = FALSE)
}

msg("Main model outputs saved.")
print(main_summary)

msg("=== Compare main model vs sensitivity analyses ===")

main_curve_cmp <- main_curve_yearly
names(main_curve_cmp)[2] <- "MainValue"

main_growth_cmp <- main_growth_yearly
names(main_growth_cmp)[2] <- "MainValue"

direct_results <- list()

for (nm in names(SensitivityFiles)) {
  msg(paste("Processing:", nm))
  
  sens_curve <- read_yearly_curve_file(SensitivityFiles[[nm]]$curve, nm, type = "curve")
  sens_growth <- read_yearly_curve_file(SensitivityFiles[[nm]]$growth, nm, type = "growth")
  
  names(sens_curve)[2] <- "SensitivityValue"
  names(sens_growth)[2] <- "SensitivityValue"
  
  res_curve <- compare_two_yearly_curves(
    main_df = main_curve_cmp,
    sens_df = sens_curve,
    analysis_name = nm,
    metric_name = "Median trajectory"
  )
  
  res_growth <- compare_two_yearly_curves(
    main_df = main_growth_cmp,
    sens_df = sens_growth,
    analysis_name = nm,
    metric_name = "Growth rate"
  )
  
  direct_results[[length(direct_results) + 1]] <- res_curve
  direct_results[[length(direct_results) + 1]] <- res_growth
}

supp_table_x <- bind_rows(direct_results)

write.csv(supp_table_x,
          file.path(OutDir, "Supplementary_Table_X_direct_comparison.csv"),
          row.names = FALSE)

msg("=== Build Supplementary Table Y ===")

summary_files <- c(
  StricterQC = "outputs/sensitivity_stricter_qc/StricterQC_manuscript_summary.csv",
  Balanced   = "outputs/sensitivity_balanced_resampling/Balanced_manuscript_summary.csv",
  SplitHalf  = "outputs/sensitivity_split_half/SplitHalf_manuscript_summary.csv",
  Bootstrap  = "outputs/sensitivity_bootstrap/Bootstrap_manuscript_summary.csv",
  LOSO       = "outputs/sensitivity_loso/LOSO_manuscript_summary.csv"
)

supp_table_y_list <- list()

for (nm in names(summary_files)) {
  f <- summary_files[[nm]]
  if (file.exists(f)) {
    tmp <- read.csv(f, stringsAsFactors = FALSE)
    tmp$AnalysisName <- nm
    supp_table_y_list[[length(supp_table_y_list) + 1]] <- tmp
  } else {
    warning(paste("Summary file not found:", f))
  }
}

supp_table_y <- bind_rows(supp_table_y_list)

write.csv(supp_table_y,
          file.path(OutDir, "Supplementary_Table_Y_sensitivity_summary.csv"),
          row.names = FALSE)

manuscript_compact <- supp_table_x %>%
  mutate(
    r_text = ifelse(is.na(r), "NA", sprintf("%.4f", r)),
    p_text = ifelse(is.na(p), "NA", format.pval(p, digits = 3, eps = .Machine$double.xmin))
  ) %>%
  select(Analysis, Metric, N_AgePoints, r, p, rmse, mean_abs_diff, max_abs_diff, r_text, p_text)

write.csv(manuscript_compact,
          file.path(OutDir, "Supplementary_Table_X_ready_for_manuscript.csv"),
          row.names = FALSE)

msg("All outputs finished.")
msg(paste("Output directory:", OutDir))
print(manuscript_compact)
