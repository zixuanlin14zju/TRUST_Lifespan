
suppressPackageStartupMessages({
  library(gamlss)
  library(pracma)   # gradient()
  library(dplyr)
  library(caret)    # createFolds
  library(moments)
  library(splines)  # bs()
})

DataPath       <- "data/trust_lifespan_hc.csv"
OutDir         <- "outputs/lifespan_model" 
OutDir2        <- "outputs/lifespan_model/data" 
if (!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)
if (!dir.exists(OutDir2)) dir.create(OutDir2, recursive = TRUE)
PhenotypeName  <- "OEF"

FontFamily     <- "Arial"
TextCex        <- 1.8      # overall text scale
AxisCex        <- 1.4      # axis tick labels
TitleCex       <- 1.5      # main title
LabCex         <- 1.5      # x/y label
PointCex       <- 1.0      # point size multiplier
LineWidth      <- 4.0      # line width for key curves

DPI            <- 600
FigW_in        <- 10
FigH_in        <- 7

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

plot_labs <- function(xlab, ylab, main){
  title(main = main, cex.main = TitleCex)
  title(xlab = xlab, ylab = ylab, cex.lab = LabCex)
  par(cex.axis = AxisCex)
}

plot_labs_small <- function(xlab, ylab, main,
                      main_cex = TitleCex*0.8,
                      main_line = 1.6){
  title(main = main, cex.main = main_cex, line = main_line)
  title(xlab = xlab, ylab = ylab, cex.lab = LabCex*0.6)
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
        names.arg = 1:17,
        col = site_cols, border = "black", las = 2)
plot_labs("Site ID", "N", "Sample size by site")
close_png()

open_png("01_02_sex_distribution_by_site.png")
gender_mat <- rbind(Female = site_stats$N_Female, Male = site_stats$N_Male)
barplot(gender_mat, beside = FALSE,
        col = c("skyblue2", "lightcoral"),
        border = "black", names.arg = 1:17, las = 2)
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
        col = age_cols, border = "black", xaxt = "n",
        cex = 0.6, xlab = "", ylab = "", main = "")
axis(1, at = 1:17, labels = 1:17, las = 2, cex.axis = AxisCex*0.7)
plot_labs("Site ID", "Age (years)", "Age by site")
close_png()

open_png("01_04_oef_by_site.png")
boxplot(phenotype ~ SiteID, data = M_HC,
        col = oef_cols, border = "black", xaxt = "n",
        cex = 0.6, xlab = "", ylab = "", main = "")
axis(1, at = 1:17, labels = 1:17, las = 2, cex.axis = AxisCex*0.7)
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

saveRDS(best_model, file = file.path(OutDir2, "best_model.rds"))

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

get_centile <- function(Cmat, q_vec, q_target) {
  idx <- which(abs(q_vec - q_target) < 1e-8)
  if (length(idx) != 1) stop("Requested quantile not found: ", q_target)
  return(as.numeric(Cmat[idx, ]))
}

P01 <- get_centile(Centiles_mean, quantiles, 0.01)
P05 <- get_centile(Centiles_mean, quantiles, 0.05)
P50 <- get_centile(Centiles_mean, quantiles, 0.50)
P95 <- get_centile(Centiles_mean, quantiles, 0.95)
P99 <- get_centile(Centiles_mean, quantiles, 0.99)

Width_90 <- P95 - P05   # 5th–95th percentile width
Width_98 <- P99 - P01   # 1st–99th percentile width

variability_df <- data.frame(
  Age = X,
  P01 = P01,
  P05 = P05,
  P50 = P50,
  P95 = P95,
  P99 = P99,
  Width_90 = Width_90,
  Width_98 = Width_98
)

variability_df$AgePeriod <- cut(
  variability_df$Age,
  breaks = c(-Inf, 18, 60, Inf),
  labels = c("Childhood/adolescence", "Adulthood", "Older age"),
  right = FALSE
)

variability_summary <- variability_df %>%
  group_by(AgePeriod) %>%
  summarise(
    Age_min = min(Age, na.rm = TRUE),
    Age_max = max(Age, na.rm = TRUE),
    
    Width90_mean = mean(Width_90, na.rm = TRUE),
    Width90_sd   = sd(Width_90, na.rm = TRUE),
    Width90_min  = min(Width_90, na.rm = TRUE),
    Width90_max  = max(Width_90, na.rm = TRUE),
    Width90_range = Width90_max - Width90_min,
    Width90_percent_change =
      100 * (Width90_max - Width90_min) / Width90_mean,
    
    Width98_mean = mean(Width_98, na.rm = TRUE),
    Width98_sd   = sd(Width_98, na.rm = TRUE),
    Width98_min  = min(Width_98, na.rm = TRUE),
    Width98_max  = max(Width_98, na.rm = TRUE),
    Width98_range = Width98_max - Width98_min,
    Width98_percent_change =
      100 * (Width98_max - Width98_min) / Width98_mean,
    
    .groups = "drop"
  )

cat("\n========== Model-derived inter-individual variability =========\n")
print(variability_summary)

write.csv(
  variability_df,
  file = file.path(OutDir2, "model_derived_centile_width_by_age.csv"),
  row.names = FALSE
)

write.csv(
  variability_summary,
  file = file.path(OutDir2, "model_derived_centile_width_summary.csv"),
  row.names = FALSE
)

open_png("02b_inter_individual_variability_centile_width.png", h = 6)

plot(
  variability_df$Age,
  variability_df$Width_90,
  type = "l",
  lwd = LineWidth,
  col = "black",
  xlab = "",
  ylab = "",
  main = ""
)

plot_labs_small(
  "Age (years)",
  paste0("P95 - P5 width of ", PhenotypeName),
  paste0("Model-derived inter-individual variability of ", PhenotypeName)
)

abline(v = c(18, 60), lty = 2)

text(
  x = c(9, 39, 76),
  y = max(variability_df$Width_90, na.rm = TRUE),
  labels = c("Childhood", "Adulthood", "Older age"),
  cex = AxisCex * 0.8,
  pos = 3
)

close_png()

make_width_df <- function(Cmat, sex_label) {
  p01 <- get_centile(Cmat, quantiles, 0.01)
  p05 <- get_centile(Cmat, quantiles, 0.05)
  p95 <- get_centile(Cmat, quantiles, 0.95)
  p99 <- get_centile(Cmat, quantiles, 0.99)
  
  data.frame(
    Age = X,
    Sex = sex_label,
    Width_90 = p95 - p05,
    Width_98 = p99 - p01
  )
}

variability_sex_df <- rbind(
  make_width_df(Centiles_female, "Female"),
  make_width_df(Centiles_male, "Male")
)

variability_sex_summary <- variability_sex_df %>%
  mutate(
    AgePeriod = cut(
      Age,
      breaks = c(-Inf, 18, 60, Inf),
      labels = c("Childhood/adolescence", "Adulthood", "Older age"),
      right = FALSE
    )
  ) %>%
  group_by(Sex, AgePeriod) %>%
  summarise(
    Width90_mean = mean(Width_90, na.rm = TRUE),
    Width90_sd   = sd(Width_90, na.rm = TRUE),
    Width90_min  = min(Width_90, na.rm = TRUE),
    Width90_max  = max(Width_90, na.rm = TRUE),
    Width90_percent_change =
      100 * (Width90_max - Width90_min) / Width90_mean,
    .groups = "drop"
  )

cat("\n========== Sex-specific inter-individual variability =========\n")
print(variability_sex_summary)

write.csv(
  variability_sex_df,
  file = file.path(OutDir2, "sex_specific_centile_width_by_age.csv"),
  row.names = FALSE
)

write.csv(
  variability_sex_summary,
  file = file.path(OutDir2, "sex_specific_centile_width_summary.csv"),
  row.names = FALSE
)

open_png("02c_sex_specific_inter_individual_variability.png", h = 6)

plot(
  variability_sex_df$Age[variability_sex_df$Sex == "Female"],
  variability_sex_df$Width_90[variability_sex_df$Sex == "Female"],
  type = "l",
  lwd = LineWidth,
  col = "indianred2",
  xlab = "",
  ylab = "",
  main = "",
  ylim = range(variability_sex_df$Width_90, na.rm = TRUE)
)

lines(
  variability_sex_df$Age[variability_sex_df$Sex == "Male"],
  variability_sex_df$Width_90[variability_sex_df$Sex == "Male"],
  lwd = LineWidth,
  col = "steelblue2"
)

plot_labs_small(
  "Age (years)",
  paste0("P95 - P5 width of ", PhenotypeName),
  paste0("Sex-specific inter-individual variability of ", PhenotypeName)
)

legend(
  "topright",
  legend = c("Female", "Male"),
  col = c("indianred2", "steelblue2"),
  lwd = LineWidth,
  bty = "n",
  cex = AxisCex
)

close_png()

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

legend("topleft", inset = c(0, -0.1),
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

open_png("04_growth_rate_overall.png", h = 6)
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

open_png("05_growth_rate_female_male.png", h = 6)

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
       inset = c(0, -0.1),
       legend = c("Female", "Male"),
       col = c("indianred2", "steelblue2"),
       lwd = LineWidth,
       bty = "n",
       cex = AxisCex)

close_png()

cat("\n========== 10-fold cross-validation z-scores (population-level) =========\n")
set.seed(123)
folds <- createFolds(M_HC$SiteID, k = 10, list = TRUE, returnTrain = FALSE)

cv_z <- rep(NA_real_, nrow(M_HC))
cv_models <- vector("list", length = 10)

for (k in 1:10) {
  cat(sprintf("Fold %d/10...\n", k))
  test_idx  <- folds[[k]]
  train_dat <- M_HC[-test_idx, ]
  test_dat  <- M_HC[test_idx, ]
  
  cv_model <- gamlss(
    phenotype ~ bs(Age, df = best_df_mu) * Sex + random(as.factor(SiteID)),
    sigma.fo  = ~ bs(Age, df = best_df_sigma) + Sex,
    nu.fo     = ~ 1,
    tau.fo    = ~ 1,
    family    = BCTo,
    data      = train_dat,
    control   = con
  )
  cv_models[[k]] <- cv_model
  
  pa <- tryCatch(
    predictAll(cv_model, newdata = test_dat, random = "zero"),
    error = function(e){
      cat(sprintf("  predictAll failed in fold %d: %s\n", k, e$message))
      return(NULL)
    }
  )
  if (is.null(pa)) {
    cv_z[test_idx] <- NA_real_
    next
  }
  
  cdf <- pBCTo(test_dat$phenotype, mu = pa$mu, sigma = pa$sigma, nu = pa$nu, tau = pa$tau)
  cdf <- pmin(pmax(cdf, 1e-10), 1 - 1e-10)
  cv_z[test_idx] <- qnorm(cdf)
  
  cat(sprintf("  Train N=%d | Test N=%d\n", nrow(train_dat), nrow(test_dat)))
}

M_HC$cv_z_score <- cv_z

cat("\n========== Traditional z-scores (population-level) =========\n")

pa_full <- tryCatch(
  predictAll(best_model, newdata = M_HC, random = "zero"),
  error = function(e){
    stop("predictAll failed for full model population-level z-scores: ", e$message)
  }
)

cdf_full <- pBCTo(M_HC$phenotype, mu = pa_full$mu, sigma = pa_full$sigma, nu = pa_full$nu, tau = pa_full$tau)
cdf_full <- pmin(pmax(cdf_full, 1e-10), 1 - 1e-10)
M_HC$traditional_z_score <- qnorm(cdf_full)

cat("\n========== Z-score summaries =========\n")
cat("CV z-score (pop-level): mean=", sprintf("%.4f", mean(M_HC$cv_z_score, na.rm = TRUE)),
    ", sd=", sprintf("%.4f", sd(M_HC$cv_z_score, na.rm = TRUE)), "\n", sep="")
cat("CV z-score range:", sprintf("%.4f", min(M_HC$cv_z_score, na.rm = TRUE)),
    "to", sprintf("%.4f", max(M_HC$cv_z_score, na.rm = TRUE)), "\n")

cat("Traditional z-score (pop-level): mean=", sprintf("%.4f", mean(M_HC$traditional_z_score, na.rm = TRUE)),
    ", sd=", sprintf("%.4f", sd(M_HC$traditional_z_score, na.rm = TRUE)), "\n", sep="")
cat("Traditional z-score range:", sprintf("%.4f", min(M_HC$traditional_z_score, na.rm = TRUE)),
    "to", sprintf("%.4f", max(M_HC$traditional_z_score, na.rm = TRUE)), "\n")

output_file <- paste0(OutDir2,"/data_with_cv_zscores.csv")
write.csv(M_HC, output_file, row.names = FALSE)

open_png("06_zscore_diagnostics.png", w = 12, h = 9)
par(mfrow = c(2, 2))

hist(M_HC$cv_z_score, breaks = 30, freq = FALSE, xlim = c(-4, 4),
     main = "", xlab = "", ylab = "")
plot_labs("z-score", "Density", paste0("CV z-score distribution (", PhenotypeName, ")"))
lines(density(M_HC$cv_z_score, na.rm = TRUE), lwd = LineWidth)
curve(dnorm(x, 0, 1), add = TRUE, lwd = max(1, LineWidth - 1), lty = 2)
legend("topright", legend = c("Empirical", "N(0,1)"),
       lwd = c(LineWidth, max(1, LineWidth - 1)), lty = c(1, 2),
       bty = "n", cex = AxisCex)

hist(M_HC$traditional_z_score, breaks = 30, freq = FALSE, xlim = c(-4, 4),
     main = "", xlab = "", ylab = "")
plot_labs("z-score", "Density", paste0("Traditional z-score distribution (", PhenotypeName, ")"))
lines(density(M_HC$traditional_z_score, na.rm = TRUE), lwd = LineWidth)
curve(dnorm(x, 0, 1), add = TRUE, lwd = max(1, LineWidth - 1), lty = 2)
legend("topright", legend = c("Empirical", "N(0,1)"),
       lwd = c(LineWidth, max(1, LineWidth - 1)), lty = c(1, 2),
       bty = "n", cex = AxisCex)

plot(M_HC$traditional_z_score, M_HC$cv_z_score,
     pch = 16, cex = PointCex,
     col = rgb(0.2, 0.4, 0.8, 0.5),
     xlim = c(-4, 4), ylim = c(-4, 4),
     main = "", xlab = "", ylab = "")
plot_labs("Traditional z-score", "CV z-score", "Traditional vs CV z-scores")
abline(a = 0, b = 1, lwd = LineWidth, lty = 2)
abline(h = 0, v = 0, lty = 2)

qqnorm(M_HC$cv_z_score,
       pch = 16, cex = PointCex,
       col = rgb(0.35, 0.35, 0.35, 0.55),
       main = "", xlab = "", ylab = "")
plot_labs("Theoretical quantiles", "Sample quantiles",
          paste0("QQ plot: CV z-score (", PhenotypeName, ")"))
qqline(M_HC$cv_z_score, lwd = LineWidth, lty = 2)

close_png()

open_png("07_cv_zscore_by_age.png", w = 10, h = 7)

plot(M_HC$Age, M_HC$cv_z_score,
     pch = 16, cex = PointCex,
     col = ifelse(M_HC$Sex == 1, "lightcoral", "skyblue2"),
     ylim = c(-6, 6),
     main = "", xlab = "", ylab = "",
     cex.axis = AxisCex)
plot_labs("Age (years)", "z-score", "z-score by age")
abline(h = c(-1.96, 0, 1.96), lty = c(2, 1, 2))
legend("topright", inset = c(0, -0.05),
       legend = c("Female", "Male"),
       col = c("lightcoral", "skyblue2"),
       pch = 16, bty = "n", cex = AxisCex*0.7)

close_png()

site_levels2 <- sort(unique(as.character(M_HC$SiteID)))
n_sites2 <- length(site_levels2)
site_cols2 <- colorRampPalette(c("#efedf5", "#54278f"))(n_sites2)

open_png("08_cv_zscore_by_site.png", w = 12, h = 7)

boxplot(cv_z_score ~ SiteID, data = M_HC,
        col = site_cols2, border = "black", xaxt = "n",
        cex = 0.6, xlab = "", ylab = "", main = "",
        ylim = c(-4, 4))
axis(1, at = 1:17, labels = 1:17, las = 2, cex.axis = AxisCex*0.7)
plot_labs("Site ID", "z-score", "z-score by site")
abline(h = c(-1.96, 0, 1.96), lty = c(2, 1, 2))

close_png()

sex_lab <- ifelse(M_HC$Sex == 1, "Female", "Male")

open_png("09_cv_zscore_by_sex.png", w = 8, h = 7)

boxplot(M_HC$cv_z_score ~ sex_lab,
        col = c("lightcoral", "skyblue2"),
        border = "black", las = 1,cex = 0.6,
        xlab = "", ylab = "", main = "",
        ylim = c(-4, 4))
plot_labs("Sex", "z-score", "z-score by sex")
abline(h = c(-1.96, 0, 1.96), lty = c(2, 1, 2))

close_png()
