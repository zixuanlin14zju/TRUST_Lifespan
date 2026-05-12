
suppressPackageStartupMessages({
  library(gamlss)
  library(pracma)   # gradient()
  library(dplyr)
  library(caret)    # createFolds
  library(moments)
  library(splines)  # bs()
})

DataPath       <- "data/trust_lifespan_hc_stricterQC.csv"
OutDir         <- "outputs/sensitivity_stricter_qc" 
if(!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
PhenotypeName  <- "OEF"

FontFamily     <- "Arial"
TextCex        <- 1.4      # overall text scale
AxisCex        <- 1.4      # axis tick labels
TitleCex       <- 1.4      # main title
LabCex         <- 1.4      # x/y label
PointCex       <- 1.0      # point size multiplier
LineWidth      <- 4.0      # line width for key curves

DPI            <- 600
FigW_in        <- 10
FigH_in        <- 6

open_png <- function(filename, w = FigW_in, h = FigH_in, res = DPI){
  f <- file.path(OutDir, filename)
  ok <- FALSE
  try({
    png(f, width = w, height = h, units = "in", res = res, type = "cairo")
    ok <- TRUE
  }, silent = TRUE)
  if (!ok) png(f, width = w, height = h, units = "in", res = res)
  par(
    family = FontFamily,
    cex    = TextCex,
    mar    = c(5, 5, 4, 2) + 0.2,
    mgp    = c(2.6, 0.9, 0),
    tcl    = -0.3,
    bty    = "l"
  )
  invisible(f)
}

close_png <- function(){
  try(dev.off(), silent = TRUE)
}

plot_labs <- function(xlab, ylab, main, sub = "with stricter QC"){
  title(main = main, cex.main = TitleCex)
  mtext(sub, side = 3, line = 0.2, cex = TitleCex*1.5, col = "gray40")
  title(xlab = xlab, ylab = ylab, cex.lab = LabCex)
  par(cex.axis = AxisCex)
}

M_HC <- read.csv(DataPath, stringsAsFactors = FALSE)

stopifnot(all(c("Age", "Sex", "SiteID", "phenotype") %in% names(M_HC)))

if (is.numeric(M_HC$Sex) || is.integer(M_HC$Sex)) {
  uniq <- sort(unique(M_HC$Sex[!is.na(M_HC$Sex)]))
  if (all(uniq %in% c(0,1))) {
    Sex_num <- as.integer(M_HC$Sex)
  } else {
    Sex_num <- as.integer(M_HC$Sex == max(uniq, na.rm = TRUE))
  }
} else {
  s <- tolower(as.character(M_HC$Sex))
  Sex_num <- ifelse(s %in% c("f", "female", "woman", "w"), 1,
                    ifelse(s %in% c("m", "male", "man"), 0, NA_integer_))
  if (mean(!is.na(Sex_num)) < 0.8) {
    s_fac <- factor(M_HC$Sex)
    levs <- levels(s_fac)
    if (length(levs) >= 2) {
      Sex_num <- as.integer(s_fac == levs[2])
    } else {
      Sex_num <- as.integer(s_fac == levs[1]) * 0  # default all 0
    }
  }
}
M_HC$Sex <- Sex_num

if (!all(na.omit(unique(M_HC$Sex)) %in% c(0,1))) {
  warning("Sex column could not be coerced cleanly to 0/1. Inspect M_HC$Sex.")
}

site_stats <- M_HC %>%
  group_by(SiteID) %>%
  summarise(
    N = n(),
    N_Female = sum(Sex == 1, na.rm = TRUE),
    N_Male   = sum(Sex == 0, na.rm = TRUE),
    Female_Percent = round(mean(Sex == 1, na.rm = TRUE) * 100, 1),
    Age_Mean = round(mean(Age, na.rm = TRUE), 1),
    Age_SD   = round(sd(Age, na.rm = TRUE), 1),
    Age_Min  = round(min(Age, na.rm = TRUE), 1),
    Age_Max  = round(max(Age, na.rm = TRUE), 1),
    Phenotype_Mean = round(mean(phenotype, na.rm = TRUE), 3),
    Phenotype_SD   = round(sd(phenotype, na.rm = TRUE), 3),
    Phenotype_Min  = round(min(phenotype, na.rm = TRUE), 3),
    Phenotype_Max  = round(max(phenotype, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  mutate(SiteID = as.character(SiteID))

cat("========== Descriptive statistics =========\n")
cat("Total N:", nrow(M_HC), "\n")
cat("Total sites:", length(unique(M_HC$SiteID)), "\n")
print(site_stats)
cat("\nOverall summary\n")
cat("Mean N per site:", round(mean(site_stats$N), 1), "\n")
cat("N range per site:", min(site_stats$N), "-", max(site_stats$N), "\n")
cat("Total females:", sum(site_stats$N_Female), "\n")
cat("Total males:", sum(site_stats$N_Male), "\n")
cat("Age mean ± SD:",
    round(mean(M_HC$Age, na.rm = TRUE), 1), "±",
    round(sd(M_HC$Age, na.rm = TRUE), 1), "years\n")
cat("Phenotype mean ± SD:",
    round(mean(M_HC$phenotype, na.rm = TRUE), 3), "±",
    round(sd(M_HC$phenotype, na.rm = TRUE), 3), "\n")

n_sites <- nrow(site_stats)
blue_palette <- colorRampPalette(c("#c6dbef", "#08306b"))(n_sites)
ord <- order(site_stats$N)          # order for mapping by sample size
site_cols <- rep(NA, n_sites)
site_cols[ord] <- blue_palette
age_cols <- colorRampPalette(c("#c7e9c0", "#00441b"))(n_sites)
oef_cols <- colorRampPalette(c("#efedf5", "#54278f"))(n_sites)

open_png("01_01_sample_size_by_site.png")
barplot(site_stats$N,
        names.arg = site_stats$SiteID,
        col = site_cols, border = "black", las = 2)
plot_labs("Site ID", "N", "Sample size by site")
close_png()

open_png("01_02_sex_distribution_by_site.png")
gender_mat <- rbind(Female = site_stats$N_Female, Male = site_stats$N_Male)
barplot(gender_mat, beside = FALSE,
        col = c("skyblue2", "lightcoral"),
        border = "black", names.arg = site_stats$SiteID, las = 2)
plot_labs("Site ID", "N", "Sex distribution by site")
legend("topright", inset = c(-0.08, 0),
       legend = c("Female", "Male"),
       fill = c("skyblue2", "lightcoral"),
       bty = "n",
       cex = AxisCex,
       xpd = NA)
close_png()

open_png("01_03_age_by_site.png")
boxplot(Age ~ SiteID, data = M_HC,
        col = age_cols, border = "black", las = 2,cex = 0.6,xlab = "", ylab = "",main = "")
plot_labs("Site ID", "Age (years)", "Age by site")
close_png()

open_png("01_04_oef_by_site.png")
boxplot(phenotype ~ SiteID, data = M_HC,
        col = oef_cols, border = "black", las = 2,cex = 0.6,xlab = "", ylab = "",main = "")
plot_labs("Site ID", PhenotypeName, paste(PhenotypeName, "by site"))
close_png()

con <- gamlss.control(n.cyc = 200)
dfs <- 2:5
best_df_sigma <- 1

models <- lapply(dfs, function(df_mu){
  gamlss(
    phenotype ~ bs(Age, df = df_mu) * Sex + random(as.factor(SiteID)),
    sigma.fo  = ~ bs(Age, df = best_df_sigma) + Sex,
    nu.fo     = ~ 1,
    tau.fo    = ~ 1,
    family    = BCTo,
    data      = M_HC,
    control   = con
  )
})

converged <- vapply(models, function(m) m$converged, logical(1))
bic_vals  <- vapply(models, function(m) m$sbc, numeric(1))
best_idx   <- which.min(bic_vals)
best_df_mu <- dfs[best_idx]
best_model <- models[[best_idx]]

cat("\n========== Model selection =========\n")
cat("Candidate df(mu):", paste(dfs, collapse = ", "), "\n")
cat("Converged:", paste(converged, collapse = ", "), "\n")
cat("BIC:", paste(round(bic_vals, 2), collapse = ", "), "\n")
cat("Best model: df(mu) =", best_df_mu, "\n")

saveRDS(best_model, file = file.path(OutDir, "best_model.rds"))

quantiles <- c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
X <- seq(-0.2, 93, length.out = 9321)

M_HC$SiteID <- as.factor(M_HC$SiteID)
site_levels_fit <- levels(M_HC$SiteID)

make_centiles <- function(sex_val, model, age_grid, q_vec, site_levels) {
  nd <- data.frame(
    Age = age_grid,
    Sex = sex_val,
    SiteID = factor(site_levels[1], levels = site_levels)
  )
  
  pa <- predictAll(model, newdata = nd, random = "zero")
  
  out <- sapply(q_vec, function(q) {
    qBCTo(q, mu = pa$mu, sigma = pa$sigma, nu = pa$nu, tau = pa$tau)
  })
  
  t(out)
}

Centiles_female <- make_centiles(
  sex_val = 1,
  model = best_model,
  age_grid = X,
  q_vec = quantiles,
  site_levels = site_levels_fit
)

Centiles_male <- make_centiles(
  sex_val = 0,
  model = best_model,
  age_grid = X,
  q_vec = quantiles,
  site_levels = site_levels_fit
)

Centiles_mean <- (Centiles_female + Centiles_male) / 2

open_png("02_centiles_overall.png")

plot(M_HC$Age, M_HC$phenotype,
     pch = 16, cex = PointCex,
     col = rgb(0.7, 0.7, 0.7, 0.5),
     xlab = "", ylab = "", main = "")

plot_labs("Age (years)", PhenotypeName,
          paste0("Lifespan centile curves of ", PhenotypeName))

for (i in seq_along(quantiles)) {
  lty_i <- if (quantiles[i] == 0.5) 1 else 2
  lwd_i <- if (quantiles[i] == 0.5) LineWidth else max(1, LineWidth - 1)
  lines(X, Centiles_mean[i, ], lwd = lwd_i, lty = lty_i)
}

close_png()

open_png("03_centiles_female_male.png")

sex_cols <- ifelse(M_HC$Sex == 1,
                   rgb(240/255, 128/255, 128/255, 0.6),  # Female
                   rgb(135/255, 206/255, 235/255, 0.6))  # Male

plot(M_HC$Age, M_HC$phenotype,
     pch = 16, cex = PointCex,
     col = sex_cols,
     xlab = "", ylab = "", main = "")

plot_labs("Age (years)", PhenotypeName,
          paste0("Lifespan centile curves of ", PhenotypeName))

idx_show <- which(quantiles %in% c(0.05, 0.5, 0.95))

for (i in idx_show) {
  lty_i <- if (quantiles[i] == 0.5) 1 else 2
  lwd_i <- if (quantiles[i] == 0.5) LineWidth else max(1, LineWidth - 1)
  
  lines(X, Centiles_female[i, ], lwd = lwd_i, lty = lty_i, col = "indianred2")
  lines(X, Centiles_male[i, ],   lwd = lwd_i, lty = lty_i, col = "steelblue2")
}

legend("topleft", inset = c(0, 0),
       legend = c("Female", "Male"),
       col = c("indianred2", "steelblue2"),
       pch = 16,
       lwd = LineWidth,
       bty = "n",
       cex = AxisCex)

close_png()

age.grid <- seq(min(M_HC$Age, na.rm = TRUE),
                max(M_HC$Age, na.rm = TRUE),
                by = 0.1)

M_HC$SiteID <- as.factor(M_HC$SiteID)
site_levels_fit <- levels(M_HC$SiteID)

make_new_mu <- function(sex_val, model, age_grid, site_levels) {
  nd <- data.frame(
    Age    = age_grid,
    Sex    = sex_val,
    SiteID = factor(site_levels[1], levels = site_levels)
  )
  
  pa <- predictAll(model, newdata = nd, random = "zero")
  return(pa$mu)
}

pred.mu_f <- make_new_mu(1, best_model, age.grid, site_levels_fit)
pred.mu_m <- make_new_mu(0, best_model, age.grid, site_levels_fit)

pred.mu_all <- (pred.mu_f + pred.mu_m) / 2

dy_dx_f   <- gradient(pred.mu_f, age.grid)
dy_dx_m   <- gradient(pred.mu_m, age.grid)
dy_dx_all <- gradient(pred.mu_all, age.grid)

y_all <- c(dy_dx_all, dy_dx_f, dy_dx_m)

y_lim_max <- max(abs(y_all), na.rm = TRUE)
y_lim <- c(0, y_lim_max)
y_ticks <- pretty(y_lim, n = 5)

open_png("04_growth_rate_overall.png", h = 6, w = 10)
plot(age.grid, dy_dx_all,
     type = "l",
     lwd = LineWidth,
     col = "black",
     ylim = y_lim,
     xlab = "", ylab = "", main = "",
     yaxt = "n")

axis(2, at = y_ticks,
     labels = sprintf("%.4f", y_ticks),
     cex.axis = AxisCex)
plot_labs("Age (years)",
          paste0("d", PhenotypeName, " / dAge"),
          paste0("Rate of change of ", PhenotypeName, " with age"))
abline(h = 0, lty = 2)
close_png()

open_png("05_growth_rate_female_male.png", h = 6, w = 10)

plot(age.grid, dy_dx_f,
     type = "l",
     lwd = LineWidth,
     col = "indianred2",
     ylim = y_lim,
     xlab = "", ylab = "", main = "",
     yaxt = "n")

axis(2, at = y_ticks,
     labels = sprintf("%.4f", y_ticks),
     cex.axis = AxisCex)

plot_labs("Age (years)",
          paste0("d", PhenotypeName, " / dAge"),
          paste0("Sex-specific rate of change of ", PhenotypeName))

lines(age.grid, dy_dx_m,
      lwd = LineWidth,
      col = "steelblue2")

abline(h = 0, lty = 2)

legend("topright",
       inset = c(0, 0),
       legend = c("Female", "Male"),
       col = c("indianred2", "steelblue2"),
       lwd = LineWidth,
       bty = "n",
       cex = AxisCex)

close_png()

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

cat("\n========== Manuscript outputs ==========\n")

idx_med <- which(quantiles == 0.5)
if (length(idx_med) != 1) stop("Median quantile 0.5 not found.")

strict_curve_all <- Centiles_mean[idx_med, ]
strict_curve_f   <- Centiles_female[idx_med, ]
strict_curve_m   <- Centiles_male[idx_med, ]

strict_curve_yearly <- sample_yearly(X, strict_curve_all, by = 1)
strict_growth_yearly <- sample_yearly(age.grid, dy_dx_all, by = 1)

strict_curve_yearly_sex <- rbind(
  data.frame(sample_yearly(X, strict_curve_f, by = 1), Sex = "Female"),
  data.frame(sample_yearly(X, strict_curve_m, by = 1), Sex = "Male")
)

strict_growth_yearly_sex <- rbind(
  data.frame(sample_yearly(age.grid, dy_dx_f, by = 1), Sex = "Female"),
  data.frame(sample_yearly(age.grid, dy_dx_m, by = 1), Sex = "Male")
)

peak_idx_all <- which.max(dy_dx_all)
trough_idx_all <- which.min(dy_dx_all)
peak_idx_f <- which.max(dy_dx_f)
trough_idx_f <- which.min(dy_dx_f)
peak_idx_m <- which.max(dy_dx_m)
trough_idx_m <- which.min(dy_dx_m)

strict_summary <- data.frame(
  Analysis = "Stricter QC",
  N_total = nrow(M_HC),
  N_sites = length(unique(M_HC$SiteID)),
  Age_min = min(M_HC$Age, na.rm = TRUE),
  Age_max = max(M_HC$Age, na.rm = TRUE),
  Best_df_mu = best_df_mu,
  PeakAge_all = age.grid[peak_idx_all],
  PeakRate_all = dy_dx_all[peak_idx_all],
  TroughAge_all = age.grid[trough_idx_all],
  TroughRate_all = dy_dx_all[trough_idx_all],
  PeakAge_female = age.grid[peak_idx_f],
  PeakRate_female = dy_dx_f[peak_idx_f],
  TroughAge_female = age.grid[trough_idx_f],
  TroughRate_female = dy_dx_f[trough_idx_f],
  PeakAge_male = age.grid[peak_idx_m],
  PeakRate_male = dy_dx_m[peak_idx_m],
  TroughAge_male = age.grid[trough_idx_m],
  TroughRate_male = dy_dx_m[trough_idx_m]
)

write.csv(strict_summary,
          file.path(OutDir, "StricterQC_manuscript_summary.csv"),
          row.names = FALSE)

write.csv(strict_curve_yearly,
          file.path(OutDir, "StricterQC_yearly_median_curve.csv"),
          row.names = FALSE)

write.csv(strict_growth_yearly,
          file.path(OutDir, "StricterQC_yearly_growth_curve.csv"),
          row.names = FALSE)

write.csv(strict_curve_yearly_sex,
          file.path(OutDir, "StricterQC_yearly_median_curve_by_sex.csv"),
          row.names = FALSE)

write.csv(strict_growth_yearly_sex,
          file.path(OutDir, "StricterQC_yearly_growth_curve_by_sex.csv"),
          row.names = FALSE)

cat("Stricter QC manuscript outputs saved.\n")
print(strict_summary)
