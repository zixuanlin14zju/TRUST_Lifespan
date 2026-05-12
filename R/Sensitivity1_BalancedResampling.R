library(dplyr)
library(gamlss)
library(gamlss.dist)
library(pracma)
library(splines)  
library(ggplot2)
set.seed(1234)

M_HC <- read.csv("data/trust_lifespan_hc.csv")
Phenotype_name <- "OEF" 

M_HC <- M_HC %>% filter(Age <= 90)
M_HC$SiteID <- factor(M_HC$SiteID)

M_HC$Sex <- as.character(M_HC$Sex)
if (!all(na.omit(unique(M_HC$Sex)) %in% c("0", "1"))) {
  stop("Sex must be coded as 0/1, with 0 = Male and 1 = Female.")
}
M_HC$Sex <- factor(M_HC$Sex, levels = c("0", "1"))

OutDir <- "outputs/sensitivity_balanced_resampling" 
if(!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)

BaseFamily   <- "Arial"
BaseSize     <- 20
TitleSize    <- 24
SubTitleSize <- 20
AxisTitle    <- 22
AxisText     <- 22
LegendText   <- 20

LineWidth    <- 1.5
MinorWidth   <- 0.9
ZeroWidth    <- 0.7
RibbonAlpha  <- 0.28

COL_F   <- "lightcoral"
COL_M   <- "skyblue2"
COL_AVG <- "mediumpurple3"

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

age_breaks  <- seq(0, 90, by = 10)
n_resample  <- 500  # change to 1000 later
quantiles   <- c(0.025, 0.5, 0.975)
age_grid    <- seq(min(M_HC$Age), max(M_HC$Age), by = 0.1)

balanced_samples <- vector("list", n_resample)

for (i in 1:n_resample) {
  tmp_data <- M_HC
  tmp_data$TempAgeGroup <- cut(tmp_data$Age, breaks = age_breaks,
                               include.lowest = TRUE, right = FALSE)
  
  min_per_group <- tmp_data %>%
    group_by(TempAgeGroup) %>%
    summarise(n = n(), .groups = "drop") %>%
    pull(n) %>% min()
  
  balanced_samples[[i]] <- tmp_data %>%
    group_by(TempAgeGroup) %>%
    sample_n(min_per_group) %>%
    ungroup() %>%
    select(-TempAgeGroup)
}

tmp_check <- M_HC
tmp_check$TempAgeGroup <- cut(
  tmp_check$Age,
  breaks = age_breaks,
  include.lowest = TRUE,
  right = FALSE
)

group_counts <- tmp_check %>%
  group_by(TempAgeGroup) %>%
  summarise(n = n(), .groups = "drop")

min_per_group <- min(group_counts$n, na.rm = TRUE)
N_groups <- nrow(group_counts)
N_final <- min_per_group * N_groups

cat("Counts in each age group:\n")
print(group_counts)

cat("\nMin per age group =", min_per_group, "\n")
cat("Number of age groups =", N_groups, "\n")
cat("Final N in each resample =", N_final, "\n")

gamlss_results <- vector("list", n_resample)

for (i in 1:n_resample) {
  tmp_data <- balanced_samples[[i]]
  
  mod <- gamlss(
    phenotype ~ bs(Age, df=2) * Sex + random(SiteID),
    sigma.fo = ~ bs(Age, df=1) + Sex,
    nu.fo = ~ 1,
    tau.fo = ~ 1,
    family = BCTo,
    data = tmp_data,
    control = gamlss.control(n.cyc = 200)
  )
  gamlss_results[[i]] <- mod
}

curve_array <- array(NA, dim = c(n_resample, length(quantiles), length(age_grid), 2))
dimnames(curve_array)[[4]] <- c("Male", "Female")

make_newdata <- function(sex_label){
  nd <- data.frame(
    Age = age_grid,
    Sex = factor(rep(sex_label, length(age_grid)), levels = c("0", "1"))
  )
  if ("SiteID" %in% names(M_HC)) {
    nd$SiteID <- factor(rep(levels(M_HC$SiteID)[1], length(age_grid)),
                        levels = levels(M_HC$SiteID))
  }
  nd
}

for (i in 1:n_resample) {
  mod <- gamlss_results[[i]]
  
  newdata_m <- make_newdata("0")
  pred_m <- predictAll(mod, newdata = newdata_m, random = "zero")
  q_matrix_m <- sapply(quantiles, function(p) {
    qBCTo(p, mu = pred_m$mu, sigma = pred_m$sigma, nu = pred_m$nu, tau = pred_m$tau)
  })
  curve_array[i, , , "Male"] <- t(q_matrix_m)
  
  newdata_f <- make_newdata("1")
  pred_f <- predictAll(mod, newdata = newdata_f, random = "zero")
  q_matrix_f <- sapply(quantiles, function(p) {
    qBCTo(p, mu = pred_f$mu, sigma = pred_f$sigma, nu = pred_f$nu, tau = pred_f$tau)
  })
  curve_array[i, , , "Female"] <- t(q_matrix_f)
}

curve_sex_avg <- apply(curve_array, c(1, 2, 3), mean, na.rm = TRUE)  # [resample, quantile, age]

overall_lower  <- apply(curve_sex_avg[, 1, ], 2, mean, na.rm = TRUE)
overall_median <- apply(curve_sex_avg[, 2, ], 2, mean, na.rm = TRUE)
overall_upper  <- apply(curve_sex_avg[, 3, ], 2, mean, na.rm = TRUE)

q_idx_lower <- which(quantiles == 0.025) 
q_idx_med <- which(quantiles == 0.5) 
q_idx_upper <- which(quantiles == 0.975) 
if (length(q_idx_lower)!=1 || length(q_idx_med)!=1 
    || length(q_idx_upper)!=1) { 
  stop("quantiles must include 0.025, 0.5, 0.975") 
} 

female_lower  <- apply(curve_array[, q_idx_lower, , "Female", drop=FALSE], 3, mean, na.rm=TRUE)
female_median <- apply(curve_array[, q_idx_med,   , "Female", drop=FALSE], 3, mean, na.rm=TRUE)
female_upper  <- apply(curve_array[, q_idx_upper, , "Female", drop=FALSE], 3, mean, na.rm=TRUE)

male_lower    <- apply(curve_array[, q_idx_lower, , "Male", drop=FALSE], 3, mean, na.rm=TRUE)
male_median   <- apply(curve_array[, q_idx_med,   , "Male", drop=FALSE], 3, mean, na.rm=TRUE)
male_upper    <- apply(curve_array[, q_idx_upper, , "Male", drop=FALSE], 3, mean, na.rm=TRUE)

growth_array <- matrix(
  NA,
  nrow = n_resample,
  ncol = length(age_grid)
)

for (i in 1:n_resample) {
  median_curve_i <- curve_sex_avg[i, 2, ]   # quantile = 0.5
  growth_array[i, ] <- gradient(median_curve_i, age_grid)
}

growth_mean <- apply(growth_array, 2, mean, na.rm = TRUE)
growth_sd   <- apply(growth_array, 2, sd,   na.rm = TRUE)
growth_lo_196 <- growth_mean - 1.96 * growth_sd
growth_hi_196 <- growth_mean + 1.96 * growth_sd

growth_female <- matrix(NA, nrow = n_resample, ncol = length(age_grid)) 
growth_male <- matrix(NA, nrow = n_resample, ncol = length(age_grid)) 
for (i in 1:n_resample) { 
  med_f <- curve_array[i, q_idx_med, , "Female"]
  med_m <- curve_array[i, q_idx_med, , "Male"]
  growth_female[i, ] <- gradient(med_f, age_grid) 
  growth_male[i, ] <- gradient(med_m, age_grid) 
} 
gmean_f <- apply(growth_female, 2, mean, na.rm = TRUE)
gsd_f   <- apply(growth_female, 2, sd,   na.rm = TRUE)
glower_f_196 <- gmean_f - 1.96 * gsd_f
gupper_f_196 <- gmean_f + 1.96 * gsd_f
gmean_m <- apply(growth_male, 2, mean, na.rm = TRUE)
gsd_m   <- apply(growth_male, 2, sd,   na.rm = TRUE)
glower_m_196 <- gmean_m - 1.96 * gsd_m
gupper_m_196 <- gmean_m + 1.96 * gsd_m

overall_median_mat <- curve_sex_avg[, q_idx_med, ]

overall_mean <- apply(overall_median_mat, 2, mean, na.rm = TRUE)
overall_sd   <- apply(overall_median_mat, 2, sd,   na.rm = TRUE)

overall_lo_196 <- overall_mean - 1.96 * overall_sd
overall_hi_196 <- overall_mean + 1.96 * overall_sd

overall_df_196 <- data.frame(
  Age   = age_grid,
  Mean  = overall_mean,
  Lo196 = overall_lo_196,
  Hi196 = overall_hi_196
)

growth_df <- data.frame(
  Age   = age_grid,
  Mean  = growth_mean,
  Lo196 = growth_lo_196,
  Hi196 = growth_hi_196
)

female_median_mat <- curve_array[, q_idx_med, , "Female"]
male_median_mat   <- curve_array[, q_idx_med, , "Male"]

female_mean <- apply(female_median_mat, 2, mean, na.rm = TRUE)
female_sd   <- apply(female_median_mat, 2, sd,   na.rm = TRUE)
female_lo_196 <- female_mean - 1.96 * female_sd
female_hi_196 <- female_mean + 1.96 * female_sd

male_mean <- apply(male_median_mat, 2, mean, na.rm = TRUE)
male_sd   <- apply(male_median_mat, 2, sd,   na.rm = TRUE)
male_lo_196 <- male_mean - 1.96 * male_sd
male_hi_196 <- male_mean + 1.96 * male_sd

sex_df_196 <- rbind(
  data.frame(Age = age_grid, Sex = "Female", Mean = female_mean, Lo196 = female_lo_196, Hi196 = female_hi_196),
  data.frame(Age = age_grid, Sex = "Male",   Mean = male_mean,   Lo196 = male_lo_196,   Hi196 = male_hi_196)
)

sex_growth_df <- rbind(
  data.frame(Age = age_grid, Sex = "Female",
             Mean = gmean_f, Lo196 = glower_f_196, Hi196 = gupper_f_196),
  data.frame(Age = age_grid, Sex = "Male",
             Mean = gmean_m, Lo196 = glower_m_196, Hi196 = gupper_m_196)
)

ylim_01_03 <- range(
  c(overall_df_196$Lo196, overall_df_196$Hi196,
    sex_df_196$Lo196, sex_df_196$Hi196),
  na.rm = TRUE
)

ylim_02_04 <- range(
  c(growth_df$Lo196, growth_df$Hi196,
    sex_growth_df$Lo196, sex_growth_df$Hi196),
  na.rm = TRUE
)

p1 <- ggplot(overall_df_196, aes(x = Age)) +
  geom_ribbon(aes(ymin = Lo196, ymax = Hi196),
              fill = COL_AVG, alpha = RibbonAlpha) +
  geom_line(aes(y = Mean),
            linewidth = LineWidth + 0.4, color = COL_AVG) +
  scale_y_continuous(limits = ylim_01_03) +
  labs(
    title = "Balanced resampling: median curve",
    x = "Age (years)",
    y = Phenotype_name
  ) +
  theme_pub()

ggsave(file.path(OutDir, "01_balancedresampling_OEF_All.png"),
       p1, width = 10, height = 6, dpi = 320)

p2 <- ggplot(growth_df, aes(x = Age, y = Mean)) +
  geom_ribbon(aes(ymin = Lo196, ymax = Hi196),
              fill = COL_AVG, alpha = RibbonAlpha) +
  geom_line(linewidth = LineWidth + 0.4, color = COL_AVG) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = ZeroWidth) +
  scale_y_continuous(limits = ylim_02_04) +
  labs(
    title = "Balanced resampling: growth rate",
    x = "Age (years)",
    y = paste0("d", Phenotype_name, " / dAge")
  ) +
  theme_pub()

ggsave(file.path(OutDir, "02_balancedresampling_growth_All.png"),
       p2, width = 10, height = 6, dpi = 320)

p3 <- ggplot(sex_df_196, aes(x = Age, y = Mean, color = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = Lo196, ymax = Hi196), alpha = RibbonAlpha) +
  geom_line(linewidth = LineWidth + 0.2) +
  scale_color_manual(values = c("Female" = COL_F, "Male" = COL_M)) +
  scale_fill_manual(values  = c("Female" = COL_F, "Male" = COL_M)) +
  scale_y_continuous(limits = ylim_01_03) +
  guides(color = guide_legend(title = NULL),
         fill  = guide_legend(title = NULL)) +
  labs(
    title = "Balanced resampling: Sex-specific median curves",
    x = "Age (years)",
    y = Phenotype_name
  ) +
  theme_pub(legend_pos = c(0.02, 0.98)) +
  theme(
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = NA, color = NA)
  )

ggsave(file.path(OutDir, "03_balancedresampling_OEF_BySex.png"),
       p3, width = 10, height = 6, dpi = 320)

p4 <- ggplot(sex_growth_df, aes(x = Age, y = Mean, color = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = Lo196, ymax = Hi196), alpha = RibbonAlpha) +
  geom_line(linewidth = LineWidth + 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = ZeroWidth) +
  scale_color_manual(values = c("Female" = COL_F, "Male" = COL_M)) +
  scale_fill_manual(values  = c("Female" = COL_F, "Male" = COL_M)) +
  scale_y_continuous(limits = ylim_02_04) +
  guides(
    color = guide_legend(title = NULL),
    fill  = guide_legend(title = NULL)
  ) +
  labs(
    title = "Balanced resampling: Sex-specific growth rates",
    x = "Age (years)",
    y = paste0("d", Phenotype_name, " / dAge")
  ) +
  theme_pub(legend_pos = c(0.98, 0.98)) +
  theme(
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = NA, color = NA)
  )

ggsave(file.path(OutDir, "04_balancedresampling_growth_BySex.png"),
       p4, width = 10, height = 6, dpi = 320)

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

cat("\n=== Step 11: Manuscript outputs ===\n")

balanced_curve_yearly <- sample_yearly(age_grid, overall_mean, by = 1)
balanced_growth_yearly <- sample_yearly(age_grid, growth_mean, by = 1)

balanced_curve_yearly_sex <- rbind(
  data.frame(sample_yearly(age_grid, female_mean, by = 1), Sex = "Female"),
  data.frame(sample_yearly(age_grid, male_mean,   by = 1), Sex = "Male")
)

balanced_growth_yearly_sex <- rbind(
  data.frame(sample_yearly(age_grid, gmean_f, by = 1), Sex = "Female"),
  data.frame(sample_yearly(age_grid, gmean_m, by = 1), Sex = "Male")
)

balanced_curve_sd_yearly <- sample_yearly(age_grid, overall_sd, by = 1)
names(balanced_curve_sd_yearly)[2] <- "SD"

balanced_growth_sd_yearly <- sample_yearly(age_grid, growth_sd, by = 1)
names(balanced_growth_sd_yearly)[2] <- "SD"

balanced_inflections <- data.frame(
  ResampleID = integer(),
  PeakAge = numeric(),
  PeakRate = numeric(),
  TroughAge = numeric(),
  TroughRate = numeric()
)

for (i in seq_len(n_resample)) {
  gr <- growth_array[i, ]
  if (all(is.na(gr))) next
  peak_idx <- which.max(gr)
  trough_idx <- which.min(gr)
  balanced_inflections <- rbind(
    balanced_inflections,
    data.frame(
      ResampleID = i,
      PeakAge = age_grid[peak_idx],
      PeakRate = gr[peak_idx],
      TroughAge = age_grid[trough_idx],
      TroughRate = gr[trough_idx]
    )
  )
}

balanced_summary <- data.frame(
  Analysis = "Balanced resampling",
  N_total = nrow(M_HC),
  N_per_resample = N_final,
  N_resamples = n_resample,
  Age_min = min(age_grid, na.rm = TRUE),
  Age_max = max(age_grid, na.rm = TRUE),
  Mean_curve_SD_across_age = mean(overall_sd, na.rm = TRUE),
  Max_curve_SD_across_age = max(overall_sd, na.rm = TRUE),
  Mean_growth_SD_across_age = mean(growth_sd, na.rm = TRUE),
  Max_growth_SD_across_age = max(growth_sd, na.rm = TRUE),
  PeakAge_mean = mean(balanced_inflections$PeakAge, na.rm = TRUE),
  PeakAge_sd = sd(balanced_inflections$PeakAge, na.rm = TRUE),
  PeakAge_q025 = quantile(balanced_inflections$PeakAge, 0.025, na.rm = TRUE),
  PeakAge_q975 = quantile(balanced_inflections$PeakAge, 0.975, na.rm = TRUE),
  TroughAge_mean = mean(balanced_inflections$TroughAge, na.rm = TRUE),
  TroughAge_sd = sd(balanced_inflections$TroughAge, na.rm = TRUE),
  TroughAge_q025 = quantile(balanced_inflections$TroughAge, 0.025, na.rm = TRUE),
  TroughAge_q975 = quantile(balanced_inflections$TroughAge, 0.975, na.rm = TRUE)
)

write.csv(balanced_summary,
          file.path(OutDir, "Balanced_manuscript_summary.csv"),
          row.names = FALSE)

write.csv(balanced_curve_yearly,
          file.path(OutDir, "Balanced_yearly_median_curve.csv"),
          row.names = FALSE)

write.csv(balanced_growth_yearly,
          file.path(OutDir, "Balanced_yearly_growth_curve.csv"),
          row.names = FALSE)

write.csv(balanced_curve_yearly_sex,
          file.path(OutDir, "Balanced_yearly_median_curve_by_sex.csv"),
          row.names = FALSE)

write.csv(balanced_growth_yearly_sex,
          file.path(OutDir, "Balanced_yearly_growth_curve_by_sex.csv"),
          row.names = FALSE)

write.csv(balanced_inflections,
          file.path(OutDir, "Balanced_inflection_points_all_resamples.csv"),
          row.names = FALSE)

cat("Balanced manuscript outputs saved.\n")
print(balanced_summary)
