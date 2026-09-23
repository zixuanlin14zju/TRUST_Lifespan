# Step26_Sensitivity_Concordance.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 26
# Split-half uses yearly First vs Second curves.
# Other resampling/QC sensitivities are compared with the main reference model.

source("R/lib/oef_utils.R")
oef_start_step("Step26")


suppressPackageStartupMessages({
  library(dplyr)
  library(gamlss)
  library(gamlss.dist)
  library(pracma)
  library(splines)
  library(readr)
})

MainDataPath <- oef_input("trust_lifespan_hc.csv")

OutDir <- oef_output("main_model_and_sensitivity_compare")
if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)

age_max <- 90
sample_age_step <- 1
main_n_cyc <- 200

main_select_df <- TRUE
main_df_candidates <- 3:5
main_mu_df <- 3
main_sigma_df <- 3

SensitivityFiles <- list(
  StricterQC = list(
    curve = oef_output("sensitivity_stricter_qc/StricterQC_yearly_median_curve.csv"),
    growth = oef_output("sensitivity_stricter_qc/StricterQC_yearly_growth_curve.csv")
  ),
  Balanced = list(
    curve = oef_output("sensitivity_balanced_resampling/Balanced_yearly_median_curve.csv"),
    growth = oef_output("sensitivity_balanced_resampling/Balanced_yearly_growth_curve.csv")
  ),
  SplitHalf = list(
    curve = oef_output("sensitivity_split_half/SplitHalf_yearly_median_curve_first_vs_second.csv"),
    growth = oef_output("sensitivity_split_half/SplitHalf_yearly_growth_first_vs_second.csv")
  ),
  Bootstrap = list(
    curve = oef_output("sensitivity_bootstrap/Bootstrap_yearly_median_curve.csv"),
    growth = oef_output("sensitivity_bootstrap/Bootstrap_yearly_growth_curve.csv")
  ),
  LOSO = list(
    curve = oef_output("sensitivity_loso/LOSO_yearly_median_curve.csv"),
    growth = oef_output("sensitivity_loso/LOSO_yearly_growth_curve.csv")
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
  out <- predict(mod, what = what, newdata = newdata, data = ref_data,
                 type = "response")
  out <- unpack_pred(out, n_target = NULL)
  if (is.null(out) || length(out) != nrow(newdata) || any(!is.finite(out))) {
    stop("Invalid newdata predictions for parameter ", what,
         "; expected ", nrow(newdata), " finite values. No padding/truncation is applied.")
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
  
  oef_fit_gamlss(
    phenotype ~ bs(Age, df = mu_df) * Sex + random(SiteID),
    sigma.fo = ~ bs(Age, df = sigma_df) + Sex,
    nu.fo = ~ 1,
    tau.fo = ~ 1,
    family = BCTo,
    data = dat,
    control = gamlss.control(n.cyc = n_cyc, trace = trace)
  )
}

fit_bcto_select_df <- function(dat, df_candidates = 3:5, sigma_df = 3, n_cyc = 200) {
  dat <- dat %>%
    dplyr::filter(!is.na(Age), !is.na(Sex), !is.na(SiteID), !is.na(phenotype))
  dat$SiteID <- factor(dat$SiteID)
  dat$Sex <- factor(as.character(dat$Sex), levels = c("0", "1"))
  
  models <- lapply(df_candidates, function(df_mu) {
    oef_fit_gamlss(
      phenotype ~ bs(Age, df = df_mu) * Sex + random(SiteID),
      sigma.fo = ~ bs(Age, df = sigma_df) + Sex,
      nu.fo = ~ 1,
      tau.fo = ~ 1,
      family = BCTo,
      data = dat,
      control = gamlss.control(n.cyc = n_cyc, trace = FALSE)
    )
  })
  
  bic_vals <- vapply(models, function(m) if (is.null(m)) NA_real_ else m$sbc, numeric(1))
  valid <- which(vapply(models, function(m) !is.null(m) && isTRUE(m$converged), logical(1)) & is.finite(bic_vals))
  if (!length(valid)) stop("No converged model for concordance reference.")
  best_idx <- valid[which.min(bic_vals[valid])]
  
  list(
    model = models[[best_idx]],
    best_df_mu = df_candidates[best_idx],
    bic_table = data.frame(df_mu = df_candidates, BIC = bic_vals)
  )
}

# BEGIN SPLIT-HALF S4 HELPERS
oef_sh_s4_number <- function(x, label) {
  text <- trimws(as.character(x))
  missing <- is.na(x) | text %in% c("", "NA", "NaN")
  value <- suppressWarnings(as.numeric(text))
  if (any(!missing & is.na(value))) {
    stop("Non-numeric value in ", label, ".", call. = FALSE)
  }
  value
}

oef_sh_s4_read_pairs <- function(file, metric) {
  if (!file.exists(file)) {
    stop("Required yearly output not found: ", file,
         "\nRun the full Step23 first, or select the directory containing its outputs.",
         call. = FALSE)
  }
  dat <- utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Age", "Value", "Half")
  if (anyDuplicated(names(dat)) || !all(required %in% names(dat))) {
    stop(basename(file), " must have unique Age, Value and Half columns.",
         call. = FALSE)
  }
  dat$Age <- oef_sh_s4_number(dat$Age, paste0(basename(file), ": Age"))
  dat$Value <- oef_sh_s4_number(dat$Value, paste0(basename(file), ": Value"))
  labels <- c("first" = "First", "first half" = "First",
              "second" = "Second", "second half" = "Second")
  dat$Half <- unname(labels[tolower(trimws(as.character(dat$Half)))])
  if (anyNA(dat$Half)) {
    stop("Unexpected or missing Half label in ", basename(file), ".", call. = FALSE)
  }
  if (any(!is.finite(dat$Age)) || any(dat$Age != floor(dat$Age))) {
    stop("Table S4 requires finite integer ages in ", basename(file), ".", call. = FALSE)
  }
  if (anyDuplicated(dat[, c("Age", "Half"), drop = FALSE])) {
    stop("Duplicate age within a half in ", basename(file),
         "; no averaging or first-row selection is applied.", call. = FALSE)
  }
  first <- dat[dat$Half == "First", c("Age", "Value"), drop = FALSE]
  second <- dat[dat$Half == "Second", c("Age", "Value"), drop = FALSE]
  first <- first[order(first$Age), , drop = FALSE]
  second <- second[order(second$Age), , drop = FALSE]
  if (nrow(first) < 3L || !identical(first$Age, second$Age)) {
    stop("First and Second must contain the same yearly ages (at least three) in ",
         basename(file), ".", call. = FALSE)
  }
  if (any(diff(first$Age) != 1)) {
    stop("The input must use consecutive 1-year age intervals: ", basename(file),
         ". No dense-grid statistics or interpolation are substituted.", call. = FALSE)
  }
  included <- is.finite(first$Value) & is.finite(second$Value)
  if (sum(included) < 3L) {
    stop("Fewer than three finite age pairs for ", metric, ".", call. = FALSE)
  }
  if (any(!included)) {
    warning(metric, ": excluding ", sum(!included),
            " non-finite age pair(s); see the agewise audit.", call. = FALSE)
  }
  difference <- first$Value - second$Value
  data.frame(
    Metric = metric,
    Age = first$Age,
    FirstHalf_Mean = first$Value,
    SecondHalf_Mean = second$Value,
    Difference_First_minus_Second = difference,
    Squared_difference = difference^2,
    Used_for_statistics = included,
    stringsAsFactors = FALSE
  )
}

oef_sh_s4_sci <- function(x) {
  out <- rep("NA", length(x))
  ok <- is.finite(x)
  out[ok] <- sub("e([+-])0+([0-9]+)$", "e\\1\\2", sprintf("%.2e", x[ok]))
  out
}

oef_sh_s4_write_numeric_csv <- function(dat, file) {
  # Use 17 significant digits for round-trip numeric audit files, with no
  # manuscript rounding. Quoted numeric cells are read back as numbers by R.
  out <- dat
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) {
      value <- out[[nm]]
      out[[nm]] <- sprintf("%.17g", value)
      out[[nm]][is.na(value)] <- NA_character_
    }
  }
  utils::write.csv(out, file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
}

oef_step26_compare_halves <- function(pairs) {
  use <- pairs$Used_for_statistics
  x <- pairs$FirstHalf_Mean[use]
  y <- pairs$SecondHalf_Mean[use]
  metric <- unique(pairs$Metric)
  if (length(metric) != 1L || length(x) < 3L) {
    stop("A single metric and at least three finite yearly pairs are required.",
         call. = FALSE)
  }
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    stop("Pearson correlation is undefined for a constant ", metric, " curve.",
         call. = FALSE)
  }
  ct <- stats::cor.test(x, y, method = "pearson", alternative = "two.sided")
  r <- unname(ct$estimate)
  data.frame(
    Analysis = "SplitHalf",
    Metric = metric,
    Comparison = "Mean First half vs mean Second half across valid random splits",
    N_AgePoints = length(x),
    r = r,
    one_minus_r = 1 - r,
    p = ct$p.value,
    rmse = sqrt(mean((x - y)^2)),
    mean_abs_diff = mean(abs(x - y)),
    max_abs_diff = max(abs(x - y)),
    stringsAsFactors = FALSE
  )
}

# END SPLIT-HALF S4 HELPERS

read_yearly_curve_file <- function(file, analysis_name, type = c("curve", "growth")) {
  type <- match.arg(type)
  if (!file.exists(file)) stop(paste("File not found:", file))
  
  dat <- read.csv(file, stringsAsFactors = FALSE)
  
  if ("Half" %in% names(dat)) {
    stop("Half-labelled curves require the dedicated First-vs-Second branch: ",
         file, ". They are not averaged before comparison with the main model.",
         call. = FALSE)
  }
  if (all(c("Age", "Value") %in% names(dat))) {
    if (anyDuplicated(dat$Age)) stop("Duplicate age rows without an explicit Half column: ", file)
    out <- dat[, c("Age", "Value")]
    names(out)[2] <- "MainLikeValue"
    return(out)
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
    Comparison = "Main model vs sensitivity mean curve",
    N_AgePoints = ct$n,
    r = ct$r,
    one_minus_r = 1 - ct$r,
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

msg("=== Concordance: main vs sensitivity; split-half First vs Second ===")

main_curve_cmp <- main_curve_yearly
names(main_curve_cmp)[2] <- "MainValue"

main_growth_cmp <- main_growth_yearly
names(main_growth_cmp)[2] <- "MainValue"

direct_results <- list()

for (nm in names(SensitivityFiles)) {
  msg(paste("Processing:", nm))

  if (identical(nm, "SplitHalf")) {
    split_curve_pairs <- oef_sh_s4_read_pairs(
      SensitivityFiles[[nm]]$curve, "Median trajectory")
    split_growth_pairs <- oef_sh_s4_read_pairs(
      SensitivityFiles[[nm]]$growth, "Growth rate")
    if (!identical(split_curve_pairs$Age, split_growth_pairs$Age)) {
      stop("Split-half median and growth files must use the same yearly age grid.",
           call. = FALSE)
    }
    res_curve <- oef_step26_compare_halves(split_curve_pairs)
    res_growth <- oef_step26_compare_halves(split_growth_pairs)
    direct_results[[length(direct_results) + 1L]] <- res_curve
    direct_results[[length(direct_results) + 1L]] <- res_growth

    split_curve_pairs$Source_file <- basename(SensitivityFiles[[nm]]$curve)
    split_curve_pairs$Source_MD5 <- unname(tools::md5sum(SensitivityFiles[[nm]]$curve))
    split_growth_pairs$Source_file <- basename(SensitivityFiles[[nm]]$growth)
    split_growth_pairs$Source_MD5 <- unname(tools::md5sum(SensitivityFiles[[nm]]$growth))
    oef_sh_s4_write_numeric_csv(
      rbind(split_curve_pairs, split_growth_pairs),
      file.path(OutDir, "SplitHalf_Table_S4_agewise_audit.csv"))
    next
  }
  
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
  StricterQC = oef_output("sensitivity_stricter_qc/StricterQC_manuscript_summary.csv"),
  Balanced   = oef_output("sensitivity_balanced_resampling/Balanced_manuscript_summary.csv"),
  SplitHalf  = oef_output("sensitivity_split_half/SplitHalf_manuscript_summary.csv"),
  Bootstrap  = oef_output("sensitivity_bootstrap/Bootstrap_manuscript_summary.csv"),
  LOSO       = oef_output("sensitivity_loso/LOSO_manuscript_summary.csv")
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
    r_text = ifelse(is.na(r), "NA", sprintf("%.9f", r)),
    p_text = ifelse(is.na(p), "NA", format.pval(p, digits = 3, eps = .Machine$double.xmin))
  ) %>%
  select(Analysis, Metric, Comparison, N_AgePoints, r, one_minus_r, p, rmse, mean_abs_diff, max_abs_diff, r_text, p_text)

write.csv(manuscript_compact,
          file.path(OutDir, "Supplementary_Table_X_ready_for_manuscript.csv"),
          row.names = FALSE)

s4_analysis_labels <- c(
  StricterQC = "Stricter QC", Balanced = "Balanced resampling",
  SplitHalf = "Split-half", Bootstrap = "Bootstrap", LOSO = "LOSO"
)
s4_order <- order(
  match(supp_table_x$Metric, c("Median trajectory", "Growth rate")),
  match(supp_table_x$Analysis, names(s4_analysis_labels))
)
s4_full <- supp_table_x[s4_order, , drop = FALSE]
s4_ready <- data.frame(
  Metric = s4_full$Metric,
  "Sensitivity analysis" = unname(s4_analysis_labels[s4_full$Analysis]),
  "Pearson's r" = sprintf("%.6f", s4_full$r),
  "1 - r" = oef_sh_s4_sci(s4_full$one_minus_r),
  "P value" = ifelse(s4_full$p < 0.0001, "<0.0001",
                      formatC(s4_full$p, format = "g", digits = 3)),
  RMSE = ifelse(s4_full$Metric == "Growth rate",
                oef_sh_s4_sci(s4_full$rmse), sprintf("%.6f", s4_full$rmse)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
rownames(s4_full) <- rownames(s4_ready) <- NULL
oef_sh_s4_write_numeric_csv(
  s4_full,
  file.path(OutDir, "Supplementary_Table_S4_resampling_full_precision.csv"))
write.csv(
  s4_ready,
  file.path(OutDir, "Supplementary_Table_S4_resampling_ready_for_manuscript.csv"),
  row.names = FALSE, fileEncoding = "UTF-8")
writeLines(c(
  "Table S4 resampling/QC rows",
  "Stricter QC, balanced resampling, bootstrap and LOSO: main model vs sensitivity mean curve.",
  "Split-half: mean First half vs mean Second half across valid random splits.",
  "All split-half statistics use the same finite pairs on the saved 1-year grid.",
  "Growth rates are read directly from the yearly growth-rate file, not re-differentiated.",
  "The original main-reference fitting/prediction and the other four comparisons are unchanged.",
  "This is a 10-row resampling/QC extract, not the full Table S4 including Hct and Ya.",
  "Hct/Ya concordance is produced by Steps 8 and 9 and is not imported by this script."
), file.path(OutDir, "Supplementary_Table_S4_resampling_notes.txt"))
msg("=== Table S4: resampling/QC rows ===")
print(s4_ready, row.names = FALSE, right = FALSE)

msg("All outputs finished.")
msg(paste("Output directory:", OutDir))
print(manuscript_compact)
