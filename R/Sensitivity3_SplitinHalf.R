suppressPackageStartupMessages({
  library(dplyr)
  library(gamlss)
  library(gamlss.dist)
  library(pracma)
  library(splines)
  library(ggplot2)
})

set.seed(1234)

InFile <- "data/trust_lifespan_hc.csv"
OutDir <- "outputs/sensitivity_split_half"

Phenotype_name <- "OEF"
N_REPS <- 500
n_points <- 1000
age_max <- 90

if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)

BaseFamily   <- "Arial"
BaseSize     <- 20
TitleSize    <- 24
SubTitleSize <- 20
AxisTitle    <- 22
AxisText     <- 22
LegendText   <- 20
LineWidth    <- 1.5
RibbonAlpha  <- 0.28

COL_F    <- "lightcoral"
COL_M    <- "skyblue2"
COL_AVG  <- "mediumpurple3"
COL_HIST <- "mediumpurple3"

theme_pub <- function(legend_pos = "none"){
  theme_classic(base_size = BaseSize, base_family = BaseFamily) +
    theme(
      plot.title    = element_text(size = TitleSize, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = SubTitleSize, color = "gray35", hjust = 0.5),
      axis.title    = element_text(size = AxisTitle),
      axis.text     = element_text(size = AxisText, color = "gray15"),
      legend.text   = element_text(size = LegendText),
      legend.title  = element_blank(),
      legend.position = legend_pos,
      panel.grid    = element_blank()
    )
}

msg <- function(...) cat(paste0(..., "\n"))

msg("=== Step 1: Load data ===")

M_HC <- read.csv(InFile)
M_HC <- M_HC %>% filter(Age <= age_max)

M_HC$SiteID <- factor(M_HC$SiteID)

M_HC$Sex <- as.character(M_HC$Sex)
if (!all(na.omit(unique(M_HC$Sex)) %in% c("0", "1"))) {
  stop("Sex must be coded as 0/1, with 0 = Male and 1 = Female.")
}
M_HC$Sex <- factor(M_HC$Sex, levels = c("0", "1"))

age_seq <- seq(min(M_HC$Age, na.rm = TRUE), max(M_HC$Age, na.rm = TRUE), length.out = n_points)

msg(sprintf("N = %d", nrow(M_HC)))
msg(sprintf("Age range = %.2f to %.2f", min(M_HC$Age, na.rm = TRUE), max(M_HC$Age, na.rm = TRUE)))
msg(sprintf("Sites = %d", nlevels(M_HC$SiteID)))

cor_results  <- rep(NA_real_, N_REPS)
rmse_results <- rep(NA_real_, N_REPS)

curve_first_mat  <- matrix(NA_real_, nrow = N_REPS, ncol = n_points)
curve_second_mat <- matrix(NA_real_, nrow = N_REPS, ncol = n_points)

growth_first_mat  <- matrix(NA_real_, nrow = N_REPS, ncol = n_points)
growth_second_mat <- matrix(NA_real_, nrow = N_REPS, ncol = n_points)

fit_and_predict_split <- function(tmp_data, age_seq) {
  data_fit <- tmp_data %>%
    filter(!is.na(Age), !is.na(Sex), !is.na(SiteID), !is.na(phenotype))
  
  data_fit$SiteID <- factor(data_fit$SiteID)
  data_fit$Sex <- factor(as.character(data_fit$Sex), levels = c("0", "1"))
  
  if (nrow(data_fit) < 20) {
    cat("  Failed: too few rows\n")
    return(NULL)
  }
  if (length(unique(data_fit$Sex)) < 2) {
    cat("  Failed: only one sex\n")
    return(NULL)
  }
  if (length(unique(data_fit$SiteID)) < 2) {
    cat("  Failed: only one site\n")
    return(NULL)
  }
  
  assign("data_fit", data_fit, envir = .GlobalEnv)
  
  mod <- tryCatch(
    gamlss(
      phenotype ~ bs(Age, df = 2) * Sex + random(SiteID),
      sigma.fo = ~ bs(Age, df = 1) + Sex,
      nu.fo    = ~ 1,
      tau.fo   = ~ 1,
      family   = BCTo,
      data     = data_fit,
      control  = gamlss.control(n.cyc = 200)
    ),
    error = function(e) {
      cat("  GAMLSS failed:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(mod)) return(NULL)
  
  site_ref <- levels(data_fit$SiteID)[1]
  
  newdata_m <- data.frame(
    Age    = age_seq,
    Sex    = factor(rep("0", length(age_seq)), levels = c("0", "1")),
    SiteID = factor(rep(site_ref, length(age_seq)), levels = levels(data_fit$SiteID))
  )
  
  newdata_f <- data.frame(
    Age    = age_seq,
    Sex    = factor(rep("1", length(age_seq)), levels = c("0", "1")),
    SiteID = factor(rep(site_ref, length(age_seq)), levels = levels(data_fit$SiteID))
  )
  
  pred_m <- tryCatch(
    predictAll(mod, newdata = newdata_m, random = "zero"),
    error = function(e) {
      cat("  predictAll male failed:", e$message, "\n")
      return(NULL)
    }
  )
  
  pred_f <- tryCatch(
    predictAll(mod, newdata = newdata_f, random = "zero"),
    error = function(e) {
      cat("  predictAll female failed:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(pred_m) || is.null(pred_f)) return(NULL)
  
  med_m <- tryCatch(
    qBCTo(0.5, mu = pred_m$mu, sigma = pred_m$sigma, nu = pred_m$nu, tau = pred_m$tau),
    error = function(e) {
      cat("  qBCTo male failed:", e$message, "\n")
      return(NULL)
    }
  )
  
  med_f <- tryCatch(
    qBCTo(0.5, mu = pred_f$mu, sigma = pred_f$sigma, nu = pred_f$nu, tau = pred_f$tau),
    error = function(e) {
      cat("  qBCTo female failed:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(med_m) || is.null(med_f)) return(NULL)
  
  list(
    median_male = med_m,
    median_female = med_f,
    median_avg = (med_m + med_f) / 2
  )
}

msg("=== Step 2: Split-half analysis ===")

for (rep in seq_len(N_REPS)) {
  cat(sprintf("\n========== Split-half run %d/%d ==========\n", rep, N_REPS))
  
  sites <- unique(M_HC$SiteID)
  M_HC_first_half  <- data.frame()
  M_HC_second_half <- data.frame()
  
  for (site in sites) {
    site_data <- M_HC[M_HC$SiteID == site, ]
    n_site <- nrow(site_data)
    
    if (n_site < 2) next
    
    half_idx <- sample(seq_len(n_site), size = floor(n_site / 2), replace = FALSE)
    M_HC_first_half  <- rbind(M_HC_first_half,  site_data[half_idx, ])
    M_HC_second_half <- rbind(M_HC_second_half, site_data[-half_idx, ])
  }
  
  cat("Processing first half...\n")
  fit_first <- fit_and_predict_split(M_HC_first_half, age_seq)
  
  cat("Processing second half...\n")
  fit_second <- fit_and_predict_split(M_HC_second_half, age_seq)
  
  if (is.null(fit_first) || is.null(fit_second)) {
    cat("One or both halves failed. Skipping this repetition.\n")
    next
  }
  
  median_centile_first  <- fit_first$median_avg
  median_centile_second <- fit_second$median_avg
  
  curve_first_mat[rep, ]  <- median_centile_first
  curve_second_mat[rep, ] <- median_centile_second
  
  growth_first_mat[rep, ]  <- gradient(median_centile_first, age_seq)
  growth_second_mat[rep, ] <- gradient(median_centile_second, age_seq)
  
  cor_val <- tryCatch(
    cor(median_centile_first, median_centile_second, use = "complete.obs"),
    error = function(e) NA_real_
  )
  
  rmse_val <- tryCatch(
    sqrt(mean((median_centile_first - median_centile_second)^2, na.rm = TRUE)),
    error = function(e) NA_real_
  )
  
  cor_results[rep]  <- cor_val
  rmse_results[rep] <- rmse_val
  
  cat(sprintf("Results: correlation = %.4f, RMSE = %.4f\n", cor_val, rmse_val))
}

msg("\n=== Step 3: Summary statistics ===")

valid_cor  <- cor_results[!is.na(cor_results)]
valid_rmse <- rmse_results[!is.na(rmse_results)]

cat(sprintf("Valid correlations: %d/%d\n", length(valid_cor), N_REPS))
cat(sprintf("Valid RMSEs: %d/%d\n", length(valid_rmse), N_REPS))

if (length(valid_cor) > 0) {
  cor_median <- median(valid_cor)
  cor_ci <- quantile(valid_cor, probs = c(0.025, 0.975))
  cor_mean <- mean(valid_cor)
  cor_sd <- sd(valid_cor)
} else {
  cor_median <- cor_ci <- cor_mean <- cor_sd <- NA
}

if (length(valid_rmse) > 0) {
  rmse_median <- median(valid_rmse)
  rmse_ci <- quantile(valid_rmse, probs = c(0.025, 0.975))
  rmse_mean <- mean(valid_rmse)
  rmse_sd <- sd(valid_rmse)
} else {
  rmse_median <- rmse_ci <- rmse_mean <- rmse_sd <- NA
}

summary_stats <- data.frame(
  Metric   = c("Pearson Correlation", "RMSE"),
  N_Valid  = c(length(valid_cor), length(valid_rmse)),
  Median   = c(cor_median, rmse_median),
  CI_Lower = c(cor_ci[1], rmse_ci[1]),
  CI_Upper = c(cor_ci[2], rmse_ci[2]),
  Mean     = c(cor_mean, rmse_mean),
  SD       = c(cor_sd, rmse_sd)
)

print(summary_stats)

valid_curve_idx <- which(!is.na(curve_first_mat[, 1]) & !is.na(curve_second_mat[, 1]))

if (length(valid_curve_idx) > 0) {
  first_valid  <- curve_first_mat[valid_curve_idx, , drop = FALSE]
  second_valid <- curve_second_mat[valid_curve_idx, , drop = FALSE]
  
  first_mean  <- apply(first_valid,  2, mean, na.rm = TRUE)
  first_sd    <- apply(first_valid,  2, sd,   na.rm = TRUE)
  second_mean <- apply(second_valid, 2, mean, na.rm = TRUE)
  second_sd   <- apply(second_valid, 2, sd,   na.rm = TRUE)
  
  curve_compare_df <- rbind(
    data.frame(Age = age_seq, Half = "First half",
               Mean = first_mean,  Lo = first_mean - 1.96 * first_sd,  Hi = first_mean + 1.96 * first_sd),
    data.frame(Age = age_seq, Half = "Second half",
               Mean = second_mean, Lo = second_mean - 1.96 * second_sd, Hi = second_mean + 1.96 * second_sd)
  )
  
  ylim_curve <- range(curve_compare_df$Lo, curve_compare_df$Hi, na.rm = TRUE)
  
  p_curve <- ggplot(curve_compare_df, aes(Age, Mean, color = Half, fill = Half)) +
    geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = RibbonAlpha) +
    geom_line(linewidth = LineWidth + 0.2) +
    scale_color_manual(values = c("First half" = COL_AVG, "Second half" = "gray50")) +
    scale_fill_manual(values  = c("First half" = COL_AVG, "Second half" = "gray50")) +
    coord_cartesian(ylim = ylim_curve) +
    guides(color = guide_legend(title = NULL),
           fill  = guide_legend(title = NULL)) +
    labs(
      title = "Split-half: median curve",
      x = "Age (years)",
      y = Phenotype_name
    ) +
    theme_pub(legend_pos = c(0.02, 0.98)) +
    theme(
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = NA, color = NA)
    )
  
  ggsave(file.path(OutDir, "09_split_half_median_curves.png"),
         p_curve, width = 10, height = 6, dpi = 320)
}

valid_growth_idx <- which(!is.na(growth_first_mat[, 1]) & !is.na(growth_second_mat[, 1]))

if (length(valid_growth_idx) > 0) {
  growth_first_valid  <- growth_first_mat[valid_growth_idx, , drop = FALSE]
  growth_second_valid <- growth_second_mat[valid_growth_idx, , drop = FALSE]
  
  first_g_mean  <- apply(growth_first_valid,  2, mean, na.rm = TRUE)
  first_g_sd    <- apply(growth_first_valid,  2, sd,   na.rm = TRUE)
  second_g_mean <- apply(growth_second_valid, 2, mean, na.rm = TRUE)
  second_g_sd   <- apply(growth_second_valid, 2, sd,   na.rm = TRUE)
  
  growth_compare_df <- rbind(
    data.frame(Age = age_seq, Half = "First half",
               Mean = first_g_mean,  Lo = first_g_mean - 1.96 * first_g_sd,  Hi = first_g_mean + 1.96 * first_g_sd),
    data.frame(Age = age_seq, Half = "Second half",
               Mean = second_g_mean, Lo = second_g_mean - 1.96 * second_g_sd, Hi = second_g_mean + 1.96 * second_g_sd)
  )
  
  ylim_growth <- range(growth_compare_df$Lo, growth_compare_df$Hi, na.rm = TRUE)
  
  p_growth <- ggplot(growth_compare_df, aes(Age, Mean, color = Half, fill = Half)) +
    geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = RibbonAlpha) +
    geom_line(linewidth = LineWidth + 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
    scale_color_manual(values = c("First half" = COL_AVG, "Second half" = "gray50")) +
    scale_fill_manual(values  = c("First half" = COL_AVG, "Second half" = "gray50")) +
    coord_cartesian(ylim = ylim_growth) +
    guides(color = guide_legend(title = NULL),
           fill  = guide_legend(title = NULL)) +
    labs(
      title = "Split-half: growth rate",
      x = "Age (years)",
      y = paste0("d", Phenotype_name, " / dAge")
    ) +
    theme_pub(legend_pos = c(0.98, 0.98)) +
    theme(
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = NA, color = NA)
    )
  
  ggsave(file.path(OutDir, "10_split_half_growth_rates.png"),
         p_growth, width = 10, height = 6, dpi = 320)
}

msg("\n=== Step 4: Histograms ===")

if (length(valid_cor) > 0) {
  df_cor <- data.frame(Value = valid_cor)
  
  p_cor <- ggplot(df_cor, aes(x = Value)) +
    geom_histogram(bins = 15, linewidth = 0.4, alpha = 0.9, fill = COL_AVG) +
    geom_vline(xintercept = cor_median, linewidth = LineWidth) +
    geom_vline(xintercept = cor_ci[1], linetype = 2, linewidth = LineWidth) +
    geom_vline(xintercept = cor_ci[2], linetype = 2, linewidth = LineWidth) +
    labs(
      title = "Split-half Pearson correlation",
      x = "Pearson r",
      y = "Count"
    ) +
    theme_pub()
  
  ggsave(file.path(OutDir, "11_split_half_hist_correlation.png"),
         p_cor, width = 10, height = 6, dpi = 320)
}

if (length(valid_rmse) > 0) {
  df_rmse <- data.frame(Value = valid_rmse)
  
  p_rmse <- ggplot(df_rmse, aes(x = Value)) +
    geom_histogram(bins = 15, linewidth = 0.4, alpha = 0.9, fill = "skyblue3") +
    geom_vline(xintercept = rmse_median, linewidth = LineWidth) +
    geom_vline(xintercept = rmse_ci[1], linetype = 2, linewidth = LineWidth) +
    geom_vline(xintercept = rmse_ci[2], linetype = 2, linewidth = LineWidth) +
    labs(
      title = "Split-half RMSE",
      x = "RMSE",
      y = "Count"
    ) +
    theme_pub()
  
  ggsave(file.path(OutDir, "12_split_half_hist_rmse.png"),
         p_rmse, width = 10, height = 6, dpi = 320)
}

msg("\n=== Step 5: Save results ===")

results_df <- data.frame(
  Repetition  = seq_len(N_REPS),
  Correlation = cor_results,
  RMSE        = rmse_results
)

write.csv(results_df,
          file.path(OutDir, "split_half_results_detailed.csv"),
          row.names = FALSE)

write.csv(summary_stats,
          file.path(OutDir, "split_half_summary_stats.csv"),
          row.names = FALSE)

if (length(valid_curve_idx) > 0) {
  write.csv(
    data.frame(
      Age = age_seq,
      FirstHalf_Mean = first_mean,
      FirstHalf_SD   = first_sd,
      FirstHalf_CI_Lower = first_mean - 1.96 * first_sd,
      FirstHalf_CI_Upper = first_mean + 1.96 * first_sd,
      SecondHalf_Mean = second_mean,
      SecondHalf_SD   = second_sd,
      SecondHalf_CI_Lower = second_mean - 1.96 * second_sd,
      SecondHalf_CI_Upper = second_mean + 1.96 * second_sd
    ),
    file.path(OutDir, "split_half_median_curves.csv"),
    row.names = FALSE
  )
}

if (length(valid_growth_idx) > 0) {
  write.csv(
    data.frame(
      Age = age_seq,
      FirstHalf_Mean = first_g_mean,
      FirstHalf_SD   = first_g_sd,
      FirstHalf_CI_Lower = first_g_mean - 1.96 * first_g_sd,
      FirstHalf_CI_Upper = first_g_mean + 1.96 * first_g_sd,
      SecondHalf_Mean = second_g_mean,
      SecondHalf_SD   = second_g_sd,
      SecondHalf_CI_Lower = second_g_mean - 1.96 * second_g_sd,
      SecondHalf_CI_Upper = second_g_mean + 1.96 * second_g_sd
    ),
    file.path(OutDir, "split_half_growth_rates.csv"),
    row.names = FALSE
  )
}

cat("\n========================================\n")
cat("Final report\n")
cat("========================================\n")
cat(sprintf("Repetitions: %d\n", N_REPS))
cat(sprintf("Valid correlations: %d (%.1f%%)\n",
            length(valid_cor), 100 * length(valid_cor) / N_REPS))
if (length(valid_cor) > 0) {
  cat(sprintf("Correlation median: %.4f (95%% CI: %.4f - %.4f)\n",
              cor_median, cor_ci[1], cor_ci[2]))
}
cat(sprintf("Valid RMSEs: %d (%.1f%%)\n",
            length(valid_rmse), 100 * length(valid_rmse) / N_REPS))
if (length(valid_rmse) > 0) {
  cat(sprintf("RMSE median: %.4f (95%% CI: %.4f - %.4f)\n",
              rmse_median, rmse_ci[1], rmse_ci[2]))
}

sample_yearly <- function(age, y, by = 1) {
  age_year <- seq(ceiling(min(age, na.rm = TRUE)),
                  floor(max(age, na.rm = TRUE)),
                  by = by)
  y_year <- approx(x = age, y = y, xout = age_year, rule = 2)$y
  data.frame(Age = age_year, Value = y_year)
}

safe_cor_test <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
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

msg("\n=== Step 9: Manuscript outputs ===")

valid_curve_idx <- which(!is.na(curve_first_mat[, 1]) & !is.na(curve_second_mat[, 1]))
valid_growth_idx <- which(!is.na(growth_first_mat[, 1]) & !is.na(growth_second_mat[, 1]))

first_valid  <- curve_first_mat[valid_curve_idx, , drop = FALSE]
second_valid <- curve_second_mat[valid_curve_idx, , drop = FALSE]

growth_first_valid  <- growth_first_mat[valid_growth_idx, , drop = FALSE]
growth_second_valid <- growth_second_mat[valid_growth_idx, , drop = FALSE]

first_mean  <- apply(first_valid,  2, mean, na.rm = TRUE)
second_mean <- apply(second_valid, 2, mean, na.rm = TRUE)
first_sd    <- apply(first_valid,  2, sd,   na.rm = TRUE)
second_sd   <- apply(second_valid, 2, sd,   na.rm = TRUE)

growth_first_mean  <- apply(growth_first_valid,  2, mean, na.rm = TRUE)
growth_second_mean <- apply(growth_second_valid, 2, mean, na.rm = TRUE)
growth_first_sd    <- apply(growth_first_valid,  2, sd,   na.rm = TRUE)
growth_second_sd   <- apply(growth_second_valid, 2, sd,   na.rm = TRUE)

split_curve_compare <- safe_cor_test(first_mean, second_mean)
split_growth_compare <- safe_cor_test(growth_first_mean, growth_second_mean)

split_summary <- data.frame(
  Analysis = "Split-half",
  N_total = nrow(M_HC),
  N_repetitions = N_REPS,
  N_valid_curve = length(valid_curve_idx),
  N_valid_growth = length(valid_growth_idx),
  Mean_half_size = round(nrow(M_HC) / 2, 2),
  Curve_r_mean = mean(valid_cor, na.rm = TRUE),
  Curve_r_sd = sd(valid_cor, na.rm = TRUE),
  Curve_r_median = median(valid_cor, na.rm = TRUE),
  Curve_r_q025 = quantile(valid_cor, 0.025, na.rm = TRUE),
  Curve_r_q975 = quantile(valid_cor, 0.975, na.rm = TRUE),
  Curve_RMSE_mean = mean(valid_rmse, na.rm = TRUE),
  Curve_RMSE_sd = sd(valid_rmse, na.rm = TRUE),
  MeanCurve_FirstVsSecond_r = split_curve_compare$r,
  MeanCurve_FirstVsSecond_p = split_curve_compare$p,
  MeanGrowth_FirstVsSecond_r = split_growth_compare$r,
  MeanGrowth_FirstVsSecond_p = split_growth_compare$p
)

split_curve_yearly <- rbind(
  data.frame(sample_yearly(age_seq, first_mean, by = 1), Half = "First"),
  data.frame(sample_yearly(age_seq, second_mean, by = 1), Half = "Second")
)

split_growth_yearly <- rbind(
  data.frame(sample_yearly(age_seq, growth_first_mean, by = 1), Half = "First"),
  data.frame(sample_yearly(age_seq, growth_second_mean, by = 1), Half = "Second")
)

split_iteration_metrics <- data.frame(
  Repetition = seq_len(N_REPS),
  Curve_r = cor_results,
  Curve_RMSE = rmse_results
)

write.csv(split_summary,
          file.path(OutDir, "SplitHalf_manuscript_summary.csv"),
          row.names = FALSE)

write.csv(split_curve_yearly,
          file.path(OutDir, "SplitHalf_yearly_median_curve_first_vs_second.csv"),
          row.names = FALSE)

write.csv(split_growth_yearly,
          file.path(OutDir, "SplitHalf_yearly_growth_first_vs_second.csv"),
          row.names = FALSE)

write.csv(split_iteration_metrics,
          file.path(OutDir, "SplitHalf_iteration_metrics.csv"),
          row.names = FALSE)

msg("Split-half manuscript outputs saved.")
print(split_summary)
