suppressPackageStartupMessages({
  library(dplyr)
  library(gamlss)
  library(gamlss.dist)
  library(pracma)
  library(ggplot2)
  library(splines)
})

set.seed(123456)

cat("=== Step 0: Data preprocessing ===\n")

M_HC <- read.csv("data/trust_lifespan_hc.csv")
Display_name <- "OEF"
output_dir <- "outputs/sensitivity_loso"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if (!("phenotype" %in% names(M_HC))) {
  stop("Column 'phenotype' not found in the input data.")
}

M_HC$SiteID <- factor(M_HC$SiteID)

M_HC$Sex <- as.character(M_HC$Sex)
if (!all(na.omit(unique(M_HC$Sex)) %in% c("0", "1"))) {
  stop("Sex must be coded as 0/1, with 0 = Male and 1 = Female.")
}
M_HC$Sex <- factor(M_HC$Sex, levels = c("0", "1"), labels = c("Male", "Female"))

site_ids <- levels(M_HC$SiteID)
n_sites <- length(site_ids)

cat(sprintf("Total sample size: %d\n", nrow(M_HC)))
cat(sprintf("Number of sites: %d\n", n_sites))
cat(sprintf("Age range: %.2f - %.2f years\n", min(M_HC$Age, na.rm = TRUE), max(M_HC$Age, na.rm = TRUE)))
cat(sprintf("Sex counts: Female=%d, Male=%d\n",
            sum(M_HC$Sex == "Female", na.rm = TRUE),
            sum(M_HC$Sex == "Male",   na.rm = TRUE)))

site_sample_sizes <- M_HC %>%
  group_by(SiteID) %>%
  summarise(n = n(), .groups = "drop")

cat("\nSample size by site:\n")
print(site_sample_sizes)

BaseFamily   <- "Arial"
BaseSize     <- 20
TitleSize    <- 24
SubTitleSize <- 20
AxisTitle    <- 22
AxisText     <- 22
LegendText   <- 20
LineWidth    <- 1.5
MinorWidth   <- 0.6
ZeroWidth    <- 0.7
RibbonAlpha  <- 0.28

COL_F    <- "lightcoral"
COL_M    <- "skyblue2"
COL_AVG  <- "mediumpurple3"
COL_GRAY <- "gray70"

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

save_png <- function(p, filename, w = 10, h = 6, dpi = 320) {
  ggsave(file.path(output_dir, filename), p, width = w, height = h, dpi = dpi)
}

cat("\n=== Step 1: Parameter settings ===\n")

age_grid <- seq(min(M_HC$Age, na.rm = TRUE), max(M_HC$Age, na.rm = TRUE), length.out = 500)

cat(sprintf("LOSO iterations: %d (one per site)\n", n_sites))
cat(sprintf("Age grid points: %d\n", length(age_grid)))
cat(sprintf("Age range: %.1f to %.1f years\n", min(age_grid), max(age_grid)))

cat("\n=== Step 2: Helper function ===\n")

fit_and_predict_loso <- function(train_data, age_grid) {
  data_fit <- train_data %>%
    filter(!is.na(Age), !is.na(Sex), !is.na(SiteID), !is.na(phenotype))
  
  data_fit$SiteID <- factor(data_fit$SiteID)
  data_fit$Sex <- factor(as.character(data_fit$Sex), levels = c("Male", "Female"))
  
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
  
  mod <- tryCatch({
    gamlss(
      phenotype ~ bs(Age, df = 2) * Sex + random(SiteID),
      sigma.fo = ~ bs(Age, df = 1) + Sex,
      nu.fo    = ~ 1,
      tau.fo   = ~ 1,
      family   = BCTo,
      data     = data_fit,
      control  = gamlss.control(n.cyc = 200, trace = FALSE)
    )
  }, error = function(e) {
    cat(sprintf("  GAMLSS failed: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(mod)) return(NULL)
  
  site_ref <- levels(data_fit$SiteID)[1]
  
  newdata_m <- data.frame(
    Age    = age_grid,
    Sex    = factor(rep("Male", length(age_grid)), levels = c("Male", "Female")),
    SiteID = factor(rep(site_ref, length(age_grid)), levels = levels(data_fit$SiteID))
  )
  
  newdata_f <- data.frame(
    Age    = age_grid,
    Sex    = factor(rep("Female", length(age_grid)), levels = c("Male", "Female")),
    SiteID = factor(rep(site_ref, length(age_grid)), levels = levels(data_fit$SiteID))
  )
  
  pred_m <- tryCatch({
    predictAll(mod, newdata = newdata_m, random = "zero")
  }, error = function(e) {
    cat(sprintf("  predictAll male failed: %s\n", e$message))
    return(NULL)
  })
  
  pred_f <- tryCatch({
    predictAll(mod, newdata = newdata_f, random = "zero")
  }, error = function(e) {
    cat(sprintf("  predictAll female failed: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(pred_m) || is.null(pred_f)) return(NULL)
  
  med_m <- tryCatch({
    qBCTo(0.5, mu = pred_m$mu, sigma = pred_m$sigma, nu = pred_m$nu, tau = pred_m$tau)
  }, error = function(e) {
    cat(sprintf("  qBCTo male failed: %s\n", e$message))
    return(NULL)
  })
  
  med_f <- tryCatch({
    qBCTo(0.5, mu = pred_f$mu, sigma = pred_f$sigma, nu = pred_f$nu, tau = pred_f$tau)
  }, error = function(e) {
    cat(sprintf("  qBCTo female failed: %s\n", e$message))
    return(NULL)
  })
  
  if (is.null(med_m) || is.null(med_f)) return(NULL)
  if (all(is.na(med_m)) || all(is.na(med_f))) {
    cat("  Failed: predicted curves are all NA\n")
    return(NULL)
  }
  
  list(
    model = mod,
    median_male = med_m,
    median_female = med_f,
    median_avg = (med_m + med_f) / 2
  )
}

cat("\n=== Step 3: Main LOSO loop ===\n")
cat("Starting LOSO analysis...\n")

loso_results <- list()

median_curves <- matrix(NA, nrow = n_sites, ncol = length(age_grid))
median_curves_male <- matrix(NA, nrow = n_sites, ncol = length(age_grid))
median_curves_female <- matrix(NA, nrow = n_sites, ncol = length(age_grid))

growth_rates <- matrix(NA, nrow = n_sites, ncol = length(age_grid))
growth_rates_male <- matrix(NA, nrow = n_sites, ncol = length(age_grid))
growth_rates_female <- matrix(NA, nrow = n_sites, ncol = length(age_grid))

param_summary <- data.frame(
  SiteID = character(),
  Mu_Intercept = numeric(),
  Mu_Terms_N = numeric(),
  Sigma_Intercept = numeric(),
  Sigma_Terms_N = numeric(),
  Nu = numeric(),
  Tau = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_len(n_sites)) {
  site_to_remove <- site_ids[i]
  
  cat(sprintf("LOSO iteration %d/%d: leaving out site %s (n=%d)\n",
              i, n_sites, site_to_remove,
              site_sample_sizes$n[site_sample_sizes$SiteID == site_to_remove]))
  
  train_data <- M_HC %>% filter(SiteID != site_to_remove)
  test_data  <- M_HC %>% filter(SiteID == site_to_remove)
  
  cat(sprintf("  Training: %d samples (%.1f%%)\n",
              nrow(train_data), 100 * nrow(train_data) / nrow(M_HC)))
  
  fit_res <- fit_and_predict_loso(train_data, age_grid)
  
  if (is.null(fit_res)) {
    cat(sprintf("  Skipping site %s\n", site_to_remove))
    next
  }
  
  mod <- fit_res$model
  median_avg <- fit_res$median_avg
  median_male <- fit_res$median_male
  median_female <- fit_res$median_female
  
  mu_params    <- tryCatch(coefficients(mod, what = "mu"),    error = function(e) NA)
  sigma_params <- tryCatch(coefficients(mod, what = "sigma"), error = function(e) NA)
  nu_params    <- tryCatch(coefficients(mod, what = "nu"),    error = function(e) NA)
  tau_params   <- tryCatch(coefficients(mod, what = "tau"),   error = function(e) NA)
  
  param_summary <- rbind(param_summary, data.frame(
    SiteID = site_to_remove,
    Mu_Intercept    = ifelse(all(is.na(mu_params)), NA, mu_params[1]),
    Mu_Terms_N      = ifelse(all(is.na(mu_params)), NA, length(mu_params)),
    Sigma_Intercept = ifelse(all(is.na(sigma_params)), NA, sigma_params[1]),
    Sigma_Terms_N   = ifelse(all(is.na(sigma_params)), NA, length(sigma_params)),
    Nu              = ifelse(all(is.na(nu_params)), NA, nu_params[1]),
    Tau             = ifelse(all(is.na(tau_params)), NA, tau_params[1])
  ))
  
  median_curves[i, ] <- median_avg
  median_curves_male[i, ] <- median_male
  median_curves_female[i, ] <- median_female
  
  growth_avg    <- gradient(median_avg, age_grid)
  growth_male   <- gradient(median_male, age_grid)
  growth_female <- gradient(median_female, age_grid)
  
  growth_rates[i, ] <- growth_avg
  growth_rates_male[i, ] <- growth_male
  growth_rates_female[i, ] <- growth_female
  
  loso_results[[i]] <- list(
    SiteID = site_to_remove,
    Model = mod,
    MedianCurve = median_avg,
    MedianCurveMale = median_male,
    MedianCurveFemale = median_female,
    GrowthRate = growth_avg,
    GrowthRateMale = growth_male,
    GrowthRateFemale = growth_female,
    TrainSize = nrow(train_data),
    TestSize = nrow(test_data)
  )
  
  cat(sprintf("  Done: site %s\n", site_to_remove))
}

successful_iterations <- sum(!is.na(median_curves[, 1]))
cat(sprintf("\nSuccessful LOSO iterations: %d/%d (%.1f%%)\n",
            successful_iterations, n_sites,
            100 * successful_iterations / n_sites))

valid_indices <- which(!is.na(median_curves[, 1]))

if (length(valid_indices) == 0) {
  stop("No valid LOSO iterations. Please inspect model fitting and prediction messages above.")
}

median_curves_valid <- median_curves[valid_indices, , drop = FALSE]
median_curves_male_valid <- median_curves_male[valid_indices, , drop = FALSE]
median_curves_female_valid <- median_curves_female[valid_indices, , drop = FALSE]

growth_rates_valid <- growth_rates[valid_indices, , drop = FALSE]
growth_rates_male_valid <- growth_rates_male[valid_indices, , drop = FALSE]
growth_rates_female_valid <- growth_rates_female[valid_indices, , drop = FALSE]

valid_site_ids <- site_ids[valid_indices]

cat("\n=== Step 4: Summary statistics ===\n")

median_mean <- apply(median_curves_valid, 2, mean, na.rm = TRUE)
median_sd   <- apply(median_curves_valid, 2, sd,   na.rm = TRUE)
median_ci_lower <- median_mean - 1.96 * median_sd
median_ci_upper <- median_mean + 1.96 * median_sd

median_male_mean <- apply(median_curves_male_valid, 2, mean, na.rm = TRUE)
median_male_sd   <- apply(median_curves_male_valid, 2, sd,   na.rm = TRUE)
median_male_ci_lower <- median_male_mean - 1.96 * median_male_sd
median_male_ci_upper <- median_male_mean + 1.96 * median_male_sd

median_female_mean <- apply(median_curves_female_valid, 2, mean, na.rm = TRUE)
median_female_sd   <- apply(median_curves_female_valid, 2, sd,   na.rm = TRUE)
median_female_ci_lower <- median_female_mean - 1.96 * median_female_sd
median_female_ci_upper <- median_female_mean + 1.96 * median_female_sd

growth_mean <- apply(growth_rates_valid, 2, mean, na.rm = TRUE)
growth_sd   <- apply(growth_rates_valid, 2, sd,   na.rm = TRUE)
growth_ci_lower <- growth_mean - 1.96 * growth_sd
growth_ci_upper <- growth_mean + 1.96 * growth_sd

growth_male_mean <- apply(growth_rates_male_valid, 2, mean, na.rm = TRUE)
growth_male_sd   <- apply(growth_rates_male_valid, 2, sd,   na.rm = TRUE)
growth_male_ci_lower <- growth_male_mean - 1.96 * growth_male_sd
growth_male_ci_upper <- growth_male_mean + 1.96 * growth_male_sd

growth_female_mean <- apply(growth_rates_female_valid, 2, mean, na.rm = TRUE)
growth_female_sd   <- apply(growth_rates_female_valid, 2, sd,   na.rm = TRUE)
growth_female_ci_lower <- growth_female_mean - 1.96 * growth_female_sd
growth_female_ci_upper <- growth_female_mean + 1.96 * growth_female_sd

median_cv <- median_sd / median_mean * 100
growth_cv <- growth_sd / pmax(abs(growth_mean), 1e-8) * 100

overall_median_sd_mean <- mean(median_sd, na.rm = TRUE)
overall_growth_sd_mean <- mean(growth_sd, na.rm = TRUE)

max_median_diff <- apply(median_curves_valid, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
max_growth_diff <- apply(growth_rates_valid, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

cat("\n=== Step 5: Single-site influence analysis ===\n")

site_deviations <- data.frame(
  SiteID = valid_site_ids,
  MeanDeviation = numeric(length(valid_site_ids)),
  MaxDeviation = numeric(length(valid_site_ids)),
  InfluenceScore = numeric(length(valid_site_ids))
)

for (i in seq_along(valid_site_ids)) {
  site_curve <- median_curves_valid[i, ]
  mean_dev <- mean(abs(site_curve - median_mean), na.rm = TRUE)
  max_dev  <- max(abs(site_curve - median_mean), na.rm = TRUE)
  
  site_sample_size <- site_sample_sizes$n[site_sample_sizes$SiteID == valid_site_ids[i]]
  influence_score <- mean_dev * log(site_sample_size + 1)
  
  site_deviations[i, ] <- c(valid_site_ids[i], mean_dev, max_dev, influence_score)
}

site_deviations$MeanDeviation  <- as.numeric(site_deviations$MeanDeviation)
site_deviations$MaxDeviation   <- as.numeric(site_deviations$MaxDeviation)
site_deviations$InfluenceScore <- as.numeric(site_deviations$InfluenceScore)

site_deviations <- site_deviations %>%
  arrange(desc(InfluenceScore))

cat("\nMost influential sites:\n")
print(head(site_deviations, 10))

param_variation <- param_summary %>%
  filter(!is.na(Mu_Intercept)) %>%
  summarise(
    Mu_Intercept_CV = sd(Mu_Intercept, na.rm = TRUE) / mean(abs(Mu_Intercept), na.rm = TRUE) * 100,
    Sigma_Intercept_CV = sd(Sigma_Intercept, na.rm = TRUE) / mean(abs(Sigma_Intercept), na.rm = TRUE) * 100,
    Nu_CV = sd(Nu, na.rm = TRUE) / mean(abs(Nu), na.rm = TRUE) * 100,
    Tau_CV = sd(Tau, na.rm = TRUE) / mean(abs(Tau), na.rm = TRUE) * 100
  )

cat("\nParameter coefficient of variation (%):\n")
print(param_variation)

cat("\n=== Step 6: Build plotting data ===\n")

spaghetti_df <- data.frame(
  Age   = rep(age_grid, nrow(median_curves_valid)),
  Value = as.vector(t(median_curves_valid)),
  Site  = rep(valid_site_ids, each = length(age_grid))
)

curve_avg_df <- data.frame(
  Age  = age_grid,
  Mean = median_mean,
  Lo   = median_ci_lower,
  Hi   = median_ci_upper
)

growth_avg_df <- data.frame(
  Age  = age_grid,
  Mean = growth_mean,
  Lo   = growth_ci_lower,
  Hi   = growth_ci_upper
)

sex_curve_df <- rbind(
  data.frame(Age = age_grid, Mean = median_male_mean,   Lo = median_male_ci_lower,   Hi = median_male_ci_upper,   Sex = "Male"),
  data.frame(Age = age_grid, Mean = median_female_mean, Lo = median_female_ci_lower, Hi = median_female_ci_upper, Sex = "Female")
)

sex_growth_df <- rbind(
  data.frame(Age = age_grid, Mean = growth_male_mean,   Lo = growth_male_ci_lower,   Hi = growth_male_ci_upper,   Sex = "Male"),
  data.frame(Age = age_grid, Mean = growth_female_mean, Lo = growth_female_ci_lower, Hi = growth_female_ci_upper, Sex = "Female")
)

cv_data <- data.frame(
  Age = age_grid,
  Median_CV = median_cv,
  GrowthRate_CV = growth_cv
)

ylim_curve_shared <- range(c(curve_avg_df$Lo, curve_avg_df$Hi,
                             sex_curve_df$Lo, sex_curve_df$Hi), na.rm = TRUE)

ylim_growth_shared <- range(c(growth_avg_df$Lo, growth_avg_df$Hi,
                              sex_growth_df$Lo, sex_growth_df$Hi), na.rm = TRUE)

cat("\n=== Step 7: Visualization ===\n")

p1 <- ggplot() +
  geom_line(
    data = spaghetti_df,
    aes(Age, Value, group = Site),
    linewidth = 0.25, alpha = 0.08, color = COL_GRAY
  ) +
  geom_ribbon(
    data = curve_avg_df,
    aes(Age, ymin = Lo, ymax = Hi),
    fill = COL_AVG, alpha = RibbonAlpha
  ) +
  geom_line(
    data = curve_avg_df,
    aes(Age, Mean),
    linewidth = LineWidth + 0.4, color = COL_AVG
  ) +
  coord_cartesian(ylim = ylim_curve_shared) +
  labs(
    title = "LOSO: median curve",
    x = "Age (years)",
    y = Display_name
  ) +
  theme_pub()

save_png(p1, "11_loso_all_curves.png", w = 10, h = 6)

p2 <- ggplot(curve_avg_df, aes(Age, Mean)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), fill = COL_AVG, alpha = RibbonAlpha) +
  geom_line(linewidth = LineWidth + 0.4, color = COL_AVG) +
  coord_cartesian(ylim = ylim_curve_shared) +
  labs(
    title = "LOSO: median curve",
    subtitle = "Mean ± 1.96 SD across LOSO iterations",
    x = "Age (years)",
    y = Display_name
  ) +
  theme_pub()

save_png(p2, "12_loso_median_curve_ci.png", w = 10, h = 6)

p3 <- ggplot(growth_avg_df, aes(Age, Mean)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), fill = COL_AVG, alpha = RibbonAlpha) +
  geom_line(linewidth = LineWidth + 0.4, color = COL_AVG) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = ZeroWidth) +
  coord_cartesian(ylim = ylim_growth_shared) +
  labs(
    title = "LOSO: growth rate",
    x = "Age (years)",
    y = paste0("d", Display_name, " / dAge")
  ) +
  theme_pub()

save_png(p3, "13_loso_growth_rate_ci.png", w = 10, h = 6)

influence_data <- site_deviations %>%
  mutate(SiteID = factor(SiteID, levels = rev(SiteID)))

top_n <- min(15, nrow(influence_data))
influence_data_top <- influence_data[1:top_n, ]

p4 <- ggplot(influence_data_top, aes(x = SiteID, y = InfluenceScore, fill = MeanDeviation)) +
  geom_col(width = 0.85) +
  scale_fill_gradient(low = "lightblue", high = "darkred", name = NULL) +
  labs(
    title = "Site influence analysis",
    subtitle = "Top 15 sites",
    x = "Site ID",
    y = "Influence score"
  ) +
  theme_pub(legend_pos = "right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_png(p4, "14_loso_site_influence.png", w = 12, h = 6)

p5 <- ggplot(sex_curve_df, aes(Age, Mean, color = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = RibbonAlpha) +
  geom_line(linewidth = LineWidth + 0.2) +
  scale_color_manual(values = c("Female" = COL_F, "Male" = COL_M)) +
  scale_fill_manual(values  = c("Female" = COL_F, "Male" = COL_M)) +
  coord_cartesian(ylim = ylim_curve_shared) +
  guides(color = guide_legend(title = NULL),
         fill  = guide_legend(title = NULL)) +
  labs(
    title = "LOSO: Sex-specific median curves",
    x = "Age (years)",
    y = Display_name
  ) +
  theme_pub(legend_pos = c(0.02, 0.98)) +
  theme(
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = NA, color = NA)
  )

save_png(p5, "15_loso_sex_specific_curves.png", w = 10, h = 6)

p6 <- ggplot(sex_growth_df, aes(Age, Mean, color = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = RibbonAlpha) +
  geom_line(linewidth = LineWidth + 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = ZeroWidth) +
  scale_color_manual(values = c("Female" = COL_F, "Male" = COL_M)) +
  scale_fill_manual(values  = c("Female" = COL_F, "Male" = COL_M)) +
  coord_cartesian(ylim = ylim_growth_shared) +
  guides(color = guide_legend(title = NULL),
         fill  = guide_legend(title = NULL)) +
  labs(
    title = "LOSO: Sex-specific growth rates",
    x = "Age (years)",
    y = paste0("d", Display_name, " / dAge")
  ) +
  theme_pub(legend_pos = c(0.98, 0.98)) +
  theme(
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = NA, color = NA)
  )

save_png(p6, "16_loso_sex_specific_growth_rates.png", w = 10, h = 6)

p7 <- ggplot(cv_data, aes(Age, Median_CV)) +
  geom_line(color = "turquoise4", linewidth = LineWidth) +
  labs(
    title = "LOSO coefficient of variation",
    subtitle = "Median curve variability across age",
    x = "Age (years)",
    y = "CoV (%)"
  ) +
  theme_pub()

save_png(p7, "17_loso_cv_median_curve.png", w = 10, h = 6)

cat("\n=== Step 8: Save result tables ===\n")

median_curve_data <- data.frame(
  Age = age_grid,
  Mean = median_mean,
  SD = median_sd,
  CI_Lower = median_ci_lower,
  CI_Upper = median_ci_upper,
  CV = median_cv
)

write.csv(median_curve_data,
          file.path(output_dir, "LOSO_median_curve_results.csv"),
          row.names = FALSE)

growth_rate_data <- data.frame(
  Age = age_grid,
  Mean = growth_mean,
  SD = growth_sd,
  CI_Lower = growth_ci_lower,
  CI_Upper = growth_ci_upper,
  CV = growth_cv
)

write.csv(growth_rate_data,
          file.path(output_dir, "LOSO_growth_rate_results.csv"),
          row.names = FALSE)

write.csv(site_deviations,
          file.path(output_dir, "LOSO_site_influence_analysis.csv"),
          row.names = FALSE)

write.csv(param_summary,
          file.path(output_dir, "LOSO_model_parameters.csv"),
          row.names = FALSE)

gender_curve_data <- data.frame(
  Age = age_grid,
  Male_Mean = median_male_mean,
  Male_SD = median_male_sd,
  Male_CI_Lower = median_male_ci_lower,
  Male_CI_Upper = median_male_ci_upper,
  Female_Mean = median_female_mean,
  Female_SD = median_female_sd,
  Female_CI_Lower = median_female_ci_lower,
  Female_CI_Upper = median_female_ci_upper
)

write.csv(gender_curve_data,
          file.path(output_dir, "LOSO_gender_specific_curves.csv"),
          row.names = FALSE)

max_sites_to_save <- min(50, nrow(median_curves_valid))
loso_curves_save <- data.frame(
  Age = rep(age_grid, max_sites_to_save),
  OEF = as.vector(t(median_curves_valid[1:max_sites_to_save, , drop = FALSE])),
  SiteID = rep(valid_site_ids[1:max_sites_to_save], each = length(age_grid))
)

write.csv(loso_curves_save,
          file.path(output_dir, "LOSO_all_loso_curves_sample.csv"),
          row.names = FALSE)

cat("\n========================================\n")
cat("LOSO analysis completed\n")
cat("========================================\n")
cat(sprintf("Successful iterations: %d/%d\n", successful_iterations, n_sites))
cat(sprintf("Mean median SD across age: %.4f\n", overall_median_sd_mean))
cat(sprintf("Mean growth-rate SD across age: %.4f\n", overall_growth_sd_mean))

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

cat("\n=== Step 9: Manuscript outputs ===\n")

loso_curve_yearly <- sample_yearly(age_grid, median_mean, by = 1)
loso_curve_yearly$SD <- sample_yearly(age_grid, median_sd, by = 1)$Value
loso_curve_yearly$CI_Lower <- sample_yearly(age_grid, median_ci_lower, by = 1)$Value
loso_curve_yearly$CI_Upper <- sample_yearly(age_grid, median_ci_upper, by = 1)$Value
loso_curve_yearly$CV <- sample_yearly(age_grid, median_cv, by = 1)$Value

loso_growth_yearly <- sample_yearly(age_grid, growth_mean, by = 1)
loso_growth_yearly$SD <- sample_yearly(age_grid, growth_sd, by = 1)$Value
loso_growth_yearly$CI_Lower <- sample_yearly(age_grid, growth_ci_lower, by = 1)$Value
loso_growth_yearly$CI_Upper <- sample_yearly(age_grid, growth_ci_upper, by = 1)$Value
loso_growth_yearly$CV <- sample_yearly(age_grid, growth_cv, by = 1)$Value

loso_inflections <- data.frame(
  SiteID = valid_site_ids,
  PeakAge = apply(growth_rates_valid, 1, function(x) age_grid[which.max(x)]),
  PeakRate = apply(growth_rates_valid, 1, function(x) max(x, na.rm = TRUE)),
  TroughAge = apply(growth_rates_valid, 1, function(x) age_grid[which.min(x)]),
  TroughRate = apply(growth_rates_valid, 1, function(x) min(x, na.rm = TRUE))
)

loso_summary <- data.frame(
  Analysis = "LOSO",
  N_total = nrow(M_HC),
  N_sites = n_sites,
  N_valid_sites = length(valid_site_ids),
  Mean_train_size = mean(sapply(loso_results[valid_indices], function(x) x$TrainSize), na.rm = TRUE),
  Mean_test_size = mean(sapply(loso_results[valid_indices], function(x) x$TestSize), na.rm = TRUE),
  Mean_curve_SD_across_age = mean(median_sd, na.rm = TRUE),
  Max_curve_SD_across_age = max(median_sd, na.rm = TRUE),
  Mean_growth_SD_across_age = mean(growth_sd, na.rm = TRUE),
  Max_growth_SD_across_age = max(growth_sd, na.rm = TRUE),
  PeakAge_mean = mean(loso_inflections$PeakAge, na.rm = TRUE),
  PeakAge_sd = sd(loso_inflections$PeakAge, na.rm = TRUE),
  PeakAge_q025 = quantile(loso_inflections$PeakAge, 0.025, na.rm = TRUE),
  PeakAge_q975 = quantile(loso_inflections$PeakAge, 0.975, na.rm = TRUE),
  TroughAge_mean = mean(loso_inflections$TroughAge, na.rm = TRUE),
  TroughAge_sd = sd(loso_inflections$TroughAge, na.rm = TRUE),
  TroughAge_q025 = quantile(loso_inflections$TroughAge, 0.025, na.rm = TRUE),
  TroughAge_q975 = quantile(loso_inflections$TroughAge, 0.975, na.rm = TRUE),
  MostInfluentialSite = site_deviations$SiteID[1],
  MaxInfluenceScore = site_deviations$InfluenceScore[1]
)

write.csv(loso_summary,
          file.path(output_dir, "LOSO_manuscript_summary.csv"),
          row.names = FALSE)

write.csv(loso_curve_yearly,
          file.path(output_dir, "LOSO_yearly_median_curve.csv"),
          row.names = FALSE)

write.csv(loso_growth_yearly,
          file.path(output_dir, "LOSO_yearly_growth_curve.csv"),
          row.names = FALSE)

write.csv(loso_inflections,
          file.path(output_dir, "LOSO_inflection_points_by_site.csv"),
          row.names = FALSE)

cat("LOSO manuscript outputs saved.\n")
print(loso_summary)
