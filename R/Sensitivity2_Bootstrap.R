
suppressPackageStartupMessages({
  library(dplyr)
  library(gamlss)
  library(gamlss.dist)
  library(pracma)    # gradient()
  library(ggplot2)
  library(cowplot)
  library(splines)   # bs()
})

set.seed(12345)  # Reproducibility

InFile <- "data/trust_lifespan_hc.csv"
OutDir <- "outputs/sensitivity_bootstrap"

Phenotype_name <- "OEF"  # column name in the CSV

n_bootstrap  <- 500        # e.g., 1000 in papers
n_age_groups <- 10
age_grid_n   <- 1000

age_df_mu    <- 2        # df for mu spline
gamlss_cyc   <- 200

BaseFamily   <- "Arial"
BaseSize     <- 20
TitleSize    <- 24
SubTitleSize <- 20
AxisTitle    <- 22
AxisText     <- 22
LegendText   <- 20
LineWidth    <- 1.5
RibbonAlpha  <- 0.18
PointSize    <- 0.7
PointAlpha   <- 0.28
PointNMax    <- 1200    

COL_F <- "lightcoral"
COL_M <- "skyblue2"
COL_AVG <- "mediumpurple3"
COL_GRAY <- "gray70"

msg <- function(...) cat(paste0(..., "\n"))

stop_if_missing <- function(dat, cols){
  miss <- setdiff(cols, names(dat))
  if(length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))
}

theme_pub <- function(legend_pos = "none"){
  theme_classic(base_size = BaseSize, base_family = BaseFamily) +
    theme(
      plot.title    = element_text(size = TitleSize, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = SubTitleSize, color = "gray35", hjust = 0.5),
      axis.title    = element_text(size = AxisTitle),
      axis.text     = element_text(size = AxisText, color = "gray15"),
      legend.title  = element_text(size = LegendText),
      legend.text   = element_text(size = LegendText),
      legend.position = legend_pos,
      panel.grid    = element_blank()
    )
}

ci_from_boot <- function(mat, z = 1.96){
  mu <- apply(mat, 2, mean, na.rm = TRUE)
  sd <- apply(mat, 2, sd,   na.rm = TRUE)
  list(mean = mu, sd = sd, lo = mu - z * sd, hi = mu + z * sd)
}

msg("=== Step 0: Data preprocessing ===")

M_HC <- read.csv(InFile)

M_HC <- M_HC %>% filter(Age <= 90)

M_HC$SiteID <- factor(M_HC$SiteID)

M_HC$Sex <- as.character(M_HC$Sex)
if (!all(na.omit(unique(M_HC$Sex)) %in% c("0", "1"))) {
  stop("Sex must be coded as 0/1, with 0 = Male and 1 = Female.")
}
M_HC$Sex <- factor(M_HC$Sex, levels = c("0","1"), labels = c("Male","Female"))

msg(sprintf("N = %d", nrow(M_HC)))
msg(sprintf("Age range: %.2f–%.2f", min(M_HC$Age, na.rm=TRUE), max(M_HC$Age, na.rm=TRUE)))
msg(sprintf("Sex counts: Female=%d, Male=%d",
            sum(M_HC$Sex=="Female", na.rm=TRUE),
            sum(M_HC$Sex=="Male",   na.rm=TRUE)))
msg(sprintf("Sites: %d", nlevels(M_HC$SiteID)))

msg("\n=== Step 1: Parameters ===")

age_breaks <- seq(min(M_HC$Age), max(M_HC$Age), length.out = n_age_groups + 1)
age_breaks[length(age_breaks)] <- age_breaks[length(age_breaks)] + 1e-8
age_grid   <- seq(min(M_HC$Age), max(M_HC$Age), length.out = age_grid_n)

msg(sprintf("Bootstrap reps: %d", n_bootstrap))
msg(sprintf("Age groups: %d", n_age_groups))
msg(sprintf("Age grid: %.1f–%.1f (n=%d)", min(age_grid), max(age_grid), length(age_grid)))

msg("\n=== Step 2: Stratified bootstrap (AgeGroup × Sex) ===")

M_HC$AgeGroup <- cut(
  M_HC$Age,
  breaks = age_breaks,
  include.lowest = TRUE,
  right = FALSE
)

if (any(is.na(M_HC$AgeGroup))) {
  warning("Some Age values were not assigned to an AgeGroup.")
}

group_counts <- M_HC %>%
  group_by(AgeGroup, Sex) %>%
  summarise(n = n(), .groups = "drop")

msg("Counts per AgeGroup × Sex:")
print(group_counts)

bootstrap_samples <- vector("list", n_bootstrap)

for (b in seq_len(n_bootstrap)) {
  if (b %% 100 == 0) msg(sprintf("  Sampling %d/%d ...", b, n_bootstrap))
  
  bootstrap_samples[[b]] <- M_HC %>%
    group_by(AgeGroup, Sex) %>%
    sample_n(size = n(), replace = TRUE) %>%
    ungroup() %>%
    select(-AgeGroup)
}

msg(sprintf("Generated %d bootstrap samples.", length(bootstrap_samples)))

msg("\n=== Step 3: Fit GAMLSS per bootstrap sample ===")

gamlss_models <- vector("list", n_bootstrap)
con <- gamlss.control(n.cyc = gamlss_cyc, trace = FALSE)

for (b in seq_len(n_bootstrap)) {
  if (b %% 100 == 0) msg(sprintf("  Fitting %d/%d ...", b, n_bootstrap))
  
  tmp_data <- bootstrap_samples[[b]]
  
  gamlss_models[[b]] <- tryCatch(
    gamlss(
      phenotype ~ bs(Age, df = age_df_mu) * Sex + random(SiteID),
      sigma.fo = ~ bs(Age, df = 1) + Sex,
      nu.fo    = ~ 1,
      tau.fo   = ~ 1,
      family   = BCTo,
      data     = tmp_data,
      control  = con
    ),
    error = function(e){
      msg(sprintf("  Model failed at bootstrap %d: %s", b, e$message))
      NULL
    }
  )
}

successful_models <- sum(!sapply(gamlss_models, is.null))
msg(sprintf("Successful models: %d/%d (%.1f%%)",
            successful_models, n_bootstrap, 100 * successful_models / n_bootstrap))

msg("\n=== Step 4: Predict median curves ===")

median_curves_male   <- matrix(NA_real_, nrow = n_bootstrap, ncol = length(age_grid))
median_curves_female <- matrix(NA_real_, nrow = n_bootstrap, ncol = length(age_grid))
median_curves_avg    <- matrix(NA_real_, nrow = n_bootstrap, ncol = length(age_grid))

site_ref <- levels(M_HC$SiteID)[1]

newdata_male <- data.frame(
  Age    = age_grid,
  Sex    = factor("Male",   levels = levels(M_HC$Sex)),
  SiteID = factor(site_ref, levels = levels(M_HC$SiteID))
)

newdata_female <- data.frame(
  Age    = age_grid,
  Sex    = factor("Female", levels = levels(M_HC$Sex)),
  SiteID = factor(site_ref, levels = levels(M_HC$SiteID))
)

for (b in seq_len(n_bootstrap)) {
  if (b %% 100 == 0) msg(sprintf("  Predicting %d/%d ...", b, n_bootstrap))
  
  mod <- gamlss_models[[b]]
  if (is.null(mod)) next
  
  pred_m <- tryCatch(predictAll(mod, newdata = newdata_male,   random = "zero"),
                     error = function(e) NULL)
  pred_f <- tryCatch(predictAll(mod, newdata = newdata_female, random = "zero"),
                     error = function(e) NULL)
  if (is.null(pred_m) || is.null(pred_f)) next
  
  med_m <- qBCTo(0.5, mu = pred_m$mu, sigma = pred_m$sigma, nu = pred_m$nu, tau = pred_m$tau)
  med_f <- qBCTo(0.5, mu = pred_f$mu, sigma = pred_f$sigma, nu = pred_f$nu, tau = pred_f$tau)
  
  median_curves_male[b, ]   <- med_m
  median_curves_female[b, ] <- med_f
  median_curves_avg[b, ]    <- (med_m + med_f) / 2
}

valid_curves <- !is.na(median_curves_avg[, 1])
msg(sprintf("Valid bootstrap curves: %d/%d", sum(valid_curves), n_bootstrap))

median_male_valid   <- median_curves_male[valid_curves, , drop = FALSE]
median_female_valid <- median_curves_female[valid_curves, , drop = FALSE]
median_avg_valid    <- median_curves_avg[valid_curves, , drop = FALSE]

msg("\n=== Step 5: Curve confidence intervals ===")

ci_avg    <- ci_from_boot(median_avg_valid)
ci_male   <- ci_from_boot(median_male_valid)
ci_female <- ci_from_boot(median_female_valid)

msg("\n=== Step 6: Growth rate (first derivative) ===")

growth_avg    <- t(apply(median_avg_valid,    1, \(x) gradient(x, age_grid)))
growth_male   <- t(apply(median_male_valid,   1, \(x) gradient(x, age_grid)))
growth_female <- t(apply(median_female_valid, 1, \(x) gradient(x, age_grid)))

ci_g_avg    <- ci_from_boot(growth_avg)
ci_g_male   <- ci_from_boot(growth_male)
ci_g_female <- ci_from_boot(growth_female)

msg("\n=== Step 7: Inflection points (peak/trough of growth rate) ===")

inflection_points <- data.frame(
  BootstrapID = integer(),
  PeakAge     = numeric(),
  PeakRate    = numeric(),
  TroughAge   = numeric(),
  TroughRate  = numeric()
)

valid_ids <- which(valid_curves)

for (k in seq_along(valid_ids)) {
  b <- valid_ids[k]
  gr <- growth_avg[k, ]
  if (all(is.na(gr))) next
  
  peak_idx   <- which.max(gr)
  trough_idx <- which.min(gr)
  
  inflection_points <- rbind(
    inflection_points,
    data.frame(
      BootstrapID = b,
      PeakAge     = age_grid[peak_idx],
      PeakRate    = gr[peak_idx],
      TroughAge   = age_grid[trough_idx],
      TroughRate  = gr[trough_idx]
    )
  )
}

if (nrow(inflection_points) > 0) {
  msg(sprintf("Peak age: %.2f ± %.2f (95%% CI: [%.2f, %.2f])",
              mean(inflection_points$PeakAge, na.rm=TRUE),
              sd(inflection_points$PeakAge,   na.rm=TRUE),
              mean(inflection_points$PeakAge, na.rm=TRUE) - 1.96 * sd(inflection_points$PeakAge, na.rm=TRUE),
              mean(inflection_points$PeakAge, na.rm=TRUE) + 1.96 * sd(inflection_points$PeakAge, na.rm=TRUE)))
}

if (!dir.exists(OutDir)) {
  dir.create(OutDir, recursive = TRUE)
  msg(sprintf("\nCreated output directory: %s", OutDir))
}

msg("\n=== Step 9: Visualization ===")

set.seed(12345)
plot_points <- M_HC %>%
  slice_sample(n = min(PointNMax, nrow(M_HC))) %>%
  mutate(Sex = factor(Sex, levels = c("Female","Male")))

spaghetti_df <- data.frame(
  Age   = rep(age_grid, nrow(median_avg_valid)),
  Value = as.vector(t(median_avg_valid)),
  Boot  = rep(seq_len(nrow(median_avg_valid)), each = length(age_grid))
)

curve_avg_df <- data.frame(
  Age  = age_grid,
  Mean = ci_avg$mean,
  Lo   = ci_avg$lo,
  Hi   = ci_avg$hi
)

growth_avg_df <- data.frame(
  Age  = age_grid,
  Mean = ci_g_avg$mean,
  Lo   = ci_g_avg$lo,
  Hi   = ci_g_avg$hi
)

sex_curve_df <- rbind(
  data.frame(Age = age_grid, Mean = ci_male$mean,   Lo = ci_male$lo,   Hi = ci_male$hi,   Sex = "Male"),
  data.frame(Age = age_grid, Mean = ci_female$mean, Lo = ci_female$lo, Hi = ci_female$hi, Sex = "Female")
)

sex_growth_df <- rbind(
  data.frame(Age = age_grid, Mean = ci_g_male$mean,   Lo = ci_g_male$lo,   Hi = ci_g_male$hi,   Sex = "Male"),
  data.frame(Age = age_grid, Mean = ci_g_female$mean, Lo = ci_g_female$lo, Hi = ci_g_female$hi, Sex = "Female")
)

ylim_05_07 <- range(
  c(curve_avg_df$Lo, curve_avg_df$Hi,
    sex_curve_df$Lo, sex_curve_df$Hi),
  na.rm = TRUE
)

ylim_06_08 <- range(
  c(growth_avg_df$Lo, growth_avg_df$Hi,
    sex_growth_df$Lo, sex_growth_df$Hi),
  na.rm = TRUE
)

p1 <- ggplot() +
  geom_line(
    data = spaghetti_df,
    aes(Age, Value, group = Boot),
    linewidth = 0.25, alpha = 0.05, color = COL_GRAY
  ) +
  geom_ribbon(
    data = curve_avg_df,
    aes(Age, ymin = Lo, ymax = Hi),
    alpha = RibbonAlpha, fill = COL_AVG
  ) +
  geom_line(
    data = curve_avg_df,
    aes(Age, Mean),
    linewidth = LineWidth + 0.4, color = COL_AVG
  ) +
  coord_cartesian(ylim = ylim_05_07) +
  labs(
    title = "Bootstrap: median curve",
    x = "Age (years)",
    y = Phenotype_name
  ) +
  theme_pub(legend_pos = "none")

ggsave(file.path(OutDir, "05_bootstrap_median_curve.png"),
       p1, width = 10, height = 6, dpi = 320)

p2 <- ggplot(growth_avg_df, aes(Age, Mean)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = RibbonAlpha, fill = COL_AVG) +
  geom_line(linewidth = LineWidth + 0.4, color = COL_AVG) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  coord_cartesian(ylim = ylim_06_08) +
  labs(
    title = "Bootstrap: growth rate",
    x = "Age (years)",
    y = paste0("d", Phenotype_name, "/dAge")
  ) +
  theme_pub(legend_pos = "none")

ggsave(file.path(OutDir, "06_bootstrap_growth_rate.png"),
       p2, width = 10, height = 6, dpi = 320)

p3 <- ggplot(sex_curve_df, aes(Age, Mean, color = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = RibbonAlpha) +
  geom_line(linewidth = LineWidth + 0.2) +
  scale_color_manual(values = c("Female" = COL_F, "Male" = COL_M)) +
  scale_fill_manual(values  = c("Female" = COL_F, "Male" = COL_M)) +
  guides(color = guide_legend(title = NULL),
         fill  = guide_legend(title = NULL)) +
  coord_cartesian(ylim = ylim_05_07) +
  labs(
    title = "Bootstrap: Sex-specific median curves",
    x = "Age (years)",
    y = Phenotype_name
  ) +
  theme_pub(legend_pos = c(0.02, 0.98)) +
  theme(
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = NA, color = NA)
  )

ggsave(file.path(OutDir, "07_Bootstrap_gender_specific_curves.png"),
       p3, width = 10, height = 6, dpi = 320)

p4 <- ggplot(sex_growth_df, aes(Age, Mean, color = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = Lo, ymax = Hi), alpha = RibbonAlpha) +
  geom_line(linewidth = LineWidth + 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  scale_color_manual(values = c("Female" = COL_F, "Male" = COL_M)) +
  scale_fill_manual(values  = c("Female" = COL_F, "Male" = COL_M)) +
  guides(color = guide_legend(title = NULL),
         fill  = guide_legend(title = NULL)) +
  coord_cartesian(ylim = ylim_06_08) +
  labs(
    title = "Bootstrap: Sex-specific growth rates",
    x = "Age (years)",
    y = paste0("d", Phenotype_name, "/dAge")
  ) +
  theme_pub(legend_pos = c(0.98, 0.98)) +
  theme(
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = NA, color = NA)
  )

ggsave(file.path(OutDir, "08_Bootstrap_gender_specific_growth_rates.png"),
       p4, width = 10, height = 6, dpi = 320)

msg("\n=== Step 10: Save results ===")

median_curve_data <- data.frame(
  Age = age_grid,
  Median_Mean = ci_avg$mean,
  Median_SD   = ci_avg$sd,
  Median_CI_Lower = ci_avg$lo,
  Median_CI_Upper = ci_avg$hi
)
write.csv(median_curve_data,
          file.path(OutDir, "median_curve_bootstrap_results.csv"),
          row.names = FALSE)

growth_rate_data <- data.frame(
  Age = age_grid,
  GrowthRate_Mean = ci_g_avg$mean,
  GrowthRate_SD   = ci_g_avg$sd,
  GrowthRate_CI_Lower = ci_g_avg$lo,
  GrowthRate_CI_Upper = ci_g_avg$hi
)
write.csv(growth_rate_data,
          file.path(OutDir, "growth_rate_bootstrap_results.csv"),
          row.names = FALSE)

gender_curve_data <- data.frame(
  Age = age_grid,
  Male_Mean   = ci_male$mean,
  Male_SD     = ci_male$sd,
  Male_CI_Lower = ci_male$lo,
  Male_CI_Upper = ci_male$hi,
  Female_Mean   = ci_female$mean,
  Female_SD     = ci_female$sd,
  Female_CI_Lower = ci_female$lo,
  Female_CI_Upper = ci_female$hi
)
write.csv(gender_curve_data,
          file.path(OutDir, "gender_specific_curves_bootstrap.csv"),
          row.names = FALSE)

if (nrow(inflection_points) > 0) {
  write.csv(inflection_points,
            file.path(OutDir, "bootstrap_inflection_points.csv"),
            row.names = FALSE)
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

msg("\n=== Step 11: Manuscript outputs ===")

bootstrap_curve_yearly <- sample_yearly(age_grid, ci_avg$mean, by = 1)
bootstrap_curve_yearly$SD <- sample_yearly(age_grid, ci_avg$sd, by = 1)$Value
bootstrap_curve_yearly$CI_Lower <- sample_yearly(age_grid, ci_avg$lo, by = 1)$Value
bootstrap_curve_yearly$CI_Upper <- sample_yearly(age_grid, ci_avg$hi, by = 1)$Value

bootstrap_growth_yearly <- sample_yearly(age_grid, ci_g_avg$mean, by = 1)
bootstrap_growth_yearly$SD <- sample_yearly(age_grid, ci_g_avg$sd, by = 1)$Value
bootstrap_growth_yearly$CI_Lower <- sample_yearly(age_grid, ci_g_avg$lo, by = 1)$Value
bootstrap_growth_yearly$CI_Upper <- sample_yearly(age_grid, ci_g_avg$hi, by = 1)$Value

bootstrap_curve_yearly_sex <- rbind(
  data.frame(sample_yearly(age_grid, ci_female$mean, by = 1), Sex = "Female"),
  data.frame(sample_yearly(age_grid, ci_male$mean,   by = 1), Sex = "Male")
)

bootstrap_growth_yearly_sex <- rbind(
  data.frame(sample_yearly(age_grid, ci_g_female$mean, by = 1), Sex = "Female"),
  data.frame(sample_yearly(age_grid, ci_g_male$mean,   by = 1), Sex = "Male")
)

bootstrap_summary <- data.frame(
  Analysis = "Bootstrap",
  N_total = nrow(M_HC),
  N_bootstrap = n_bootstrap,
  N_valid = sum(valid_curves, na.rm = TRUE),
  Age_min = min(age_grid, na.rm = TRUE),
  Age_max = max(age_grid, na.rm = TRUE),
  Mean_curve_SD_across_age = mean(ci_avg$sd, na.rm = TRUE),
  Max_curve_SD_across_age = max(ci_avg$sd, na.rm = TRUE),
  Mean_growth_SD_across_age = mean(ci_g_avg$sd, na.rm = TRUE),
  Max_growth_SD_across_age = max(ci_g_avg$sd, na.rm = TRUE),
  PeakAge_mean = mean(inflection_points$PeakAge, na.rm = TRUE),
  PeakAge_sd = sd(inflection_points$PeakAge, na.rm = TRUE),
  PeakAge_q025 = quantile(inflection_points$PeakAge, 0.025, na.rm = TRUE),
  PeakAge_q975 = quantile(inflection_points$PeakAge, 0.975, na.rm = TRUE),
  TroughAge_mean = mean(inflection_points$TroughAge, na.rm = TRUE),
  TroughAge_sd = sd(inflection_points$TroughAge, na.rm = TRUE),
  TroughAge_q025 = quantile(inflection_points$TroughAge, 0.025, na.rm = TRUE),
  TroughAge_q975 = quantile(inflection_points$TroughAge, 0.975, na.rm = TRUE)
)

write.csv(bootstrap_summary,
          file.path(OutDir, "Bootstrap_manuscript_summary.csv"),
          row.names = FALSE)

write.csv(bootstrap_curve_yearly,
          file.path(OutDir, "Bootstrap_yearly_median_curve.csv"),
          row.names = FALSE)

write.csv(bootstrap_growth_yearly,
          file.path(OutDir, "Bootstrap_yearly_growth_curve.csv"),
          row.names = FALSE)

write.csv(bootstrap_curve_yearly_sex,
          file.path(OutDir, "Bootstrap_yearly_median_curve_by_sex.csv"),
          row.names = FALSE)

write.csv(bootstrap_growth_yearly_sex,
          file.path(OutDir, "Bootstrap_yearly_growth_curve_by_sex.csv"),
          row.names = FALSE)

write.csv(inflection_points,
          file.path(OutDir, "Bootstrap_inflection_points_all_iterations.csv"),
          row.names = FALSE)

msg("Bootstrap manuscript outputs saved.")
print(bootstrap_summary)
