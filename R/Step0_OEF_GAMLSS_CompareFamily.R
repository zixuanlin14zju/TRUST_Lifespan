
suppressPackageStartupMessages({
  library(gamlss)
  library(gamlss.dist)
  library(dplyr)
  library(splines)
  library(ggplot2)
})

DataPath <- "data/trust_lifespan_hc.csv"
OutDir2  <- "outputs/model_selection"

if (!dir.exists(OutDir2)) dir.create(OutDir2, recursive = TRUE)

set.seed(123)

dfs_mu    <- 2:6
dfs_sigma <- 1:4

family_list <- list(
  NO    = NO,
  LO    = LO,
  BCCG  = BCCG,
  BCPE  = BCPE,
  BCTo  = BCTo,
  SHASH = SHASH,
  JSU   = JSU
)

con <- gamlss.control(n.cyc = 200, trace = FALSE)

M_HC <- read.csv(DataPath, stringsAsFactors = FALSE)

stopifnot(all(c("Age", "Sex", "SiteID", "phenotype") %in% names(M_HC)))

if (is.numeric(M_HC$Sex) || is.integer(M_HC$Sex)) {
  uniq <- sort(unique(M_HC$Sex[!is.na(M_HC$Sex)]))
  if (all(uniq %in% c(0, 1))) {
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
      Sex_num <- rep(0L, nrow(M_HC))
    }
  }
}
M_HC$Sex <- Sex_num

M_HC <- M_HC %>%
  filter(
    !is.na(Age),
    !is.na(Sex),
    !is.na(SiteID),
    !is.na(phenotype)
  )

cat("N used for model comparison:", nrow(M_HC), "\n")
cat("Any phenotype <= 0 ? ", any(M_HC$phenotype <= 0, na.rm = TRUE), "\n")

fit_gamlss_safe <- function(data, df_mu, df_sigma, fam_fun, control_obj) {
  out <- tryCatch(
    gamlss(
      phenotype ~ bs(Age, df = df_mu) + Sex + random(as.factor(SiteID)),
      sigma.fo  = ~ bs(Age, df = df_sigma) + Sex,
      nu.fo     = ~ 1,
      tau.fo    = ~ 1,
      family    = fam_fun,
      data      = data,
      control   = control_obj
    ),
    error = function(e) e
  )
  return(out)
}

grid_df <- expand.grid(
  Family   = names(family_list),
  DF_Mu    = dfs_mu,
  DF_Sigma = dfs_sigma,
  stringsAsFactors = FALSE
)

results_list <- vector("list", nrow(grid_df))
model_store  <- list()

for (i in seq_len(nrow(grid_df))) {
  fam_name <- grid_df$Family[i]
  df_mu    <- grid_df$DF_Mu[i]
  df_sigma <- grid_df$DF_Sigma[i]
  fam_fun  <- family_list[[fam_name]]
  
  cat(sprintf("[%d / %d] Family=%s | df_mu=%d | df_sigma=%d\n",
              i, nrow(grid_df), fam_name, df_mu, df_sigma))
  
  fit <- fit_gamlss_safe(
    data = M_HC,
    df_mu = df_mu,
    df_sigma = df_sigma,
    fam_fun = fam_fun,
    control_obj = con
  )
  
  if (inherits(fit, "error")) {
    cat("  Model failed:", fit$message, "\n")
    
    results_list[[i]] <- data.frame(
      Family = fam_name,
      DF_Mu = df_mu,
      DF_Sigma = df_sigma,
      Converged = FALSE,
      AIC = NA_real_,
      BIC = NA_real_,
      GAIC2 = NA_real_,
      GlobalDeviance = NA_real_,
      stringsAsFactors = FALSE
    )
    next
  }
  
  conv_flag <- isTRUE(fit$converged)
  if (!conv_flag) cat("  Model did not converge.\n")
  
  aic_val  <- tryCatch(AIC(fit), error = function(e) NA_real_)
  bic_val  <- tryCatch(BIC(fit), error = function(e) NA_real_)
  gaic_val <- tryCatch(GAIC(fit, k = 2), error = function(e) NA_real_)
  gd_val   <- tryCatch(deviance(fit), error = function(e) NA_real_)
  
  results_list[[i]] <- data.frame(
    Family = fam_name,
    DF_Mu = df_mu,
    DF_Sigma = df_sigma,
    Converged = conv_flag,
    AIC = aic_val,
    BIC = bic_val,
    GAIC2 = gaic_val,
    GlobalDeviance = gd_val,
    stringsAsFactors = FALSE
  )
  
  if (conv_flag) {
    key <- paste0(fam_name, "_dfmu", df_mu, "_dfsigma", df_sigma)
    model_store[[key]] <- fit
  }
}

results_df <- bind_rows(results_list)

results_ranked_bic <- results_df %>%
  filter(Converged) %>%
  arrange(BIC, AIC)

results_ranked_aic <- results_df %>%
  filter(Converged) %>%
  arrange(AIC, BIC)

cat("\n==================== Top models by BIC ====================\n")
print(head(results_ranked_bic, 20))

cat("\n==================== Top models by AIC ====================\n")
print(head(results_ranked_aic, 20))

write.csv(
  results_df,
  file.path(OutDir2, "gamlss_family_dfmu_dfsigma_AIC_BIC_all.csv"),
  row.names = FALSE
)

write.csv(
  results_ranked_bic,
  file.path(OutDir2, "gamlss_family_dfmu_dfsigma_ranked_by_BIC.csv"),
  row.names = FALSE
)

write.csv(
  results_ranked_aic,
  file.path(OutDir2, "gamlss_family_dfmu_dfsigma_ranked_by_AIC.csv"),
  row.names = FALSE
)

if (nrow(results_ranked_bic) == 0) {
  stop("No converged model found.")
}

best_row_bic <- results_ranked_bic[1, ]
best_key_bic <- paste0(best_row_bic$Family, "_dfmu", best_row_bic$DF_Mu, "_dfsigma", best_row_bic$DF_Sigma)
best_model_bic <- model_store[[best_key_bic]]

saveRDS(
  best_model_bic,
  file = file.path(OutDir2, "best_model_by_BIC.rds")
)

write.csv(
  best_row_bic,
  file.path(OutDir2, "best_model_by_BIC_summary.csv"),
  row.names = FALSE
)

best_row_aic <- results_ranked_aic[1, ]
best_key_aic <- paste0(best_row_aic$Family, "_dfmu", best_row_aic$DF_Mu, "_dfsigma", best_row_aic$DF_Sigma)
best_model_aic <- model_store[[best_key_aic]]

saveRDS(
  best_model_aic,
  file = file.path(OutDir2, "best_model_by_AIC.rds")
)

write.csv(
  best_row_aic,
  file.path(OutDir2, "best_model_by_AIC_summary.csv"),
  row.names = FALSE
)

cat("\n==================== Best model by BIC ====================\n")
cat("Family   :", best_row_bic$Family, "\n")
cat("df_mu    :", best_row_bic$DF_Mu, "\n")
cat("df_sigma :", best_row_bic$DF_Sigma, "\n")
cat("AIC      :", round(best_row_bic$AIC, 2), "\n")
cat("BIC      :", round(best_row_bic$BIC, 2), "\n")
cat("GAIC2    :", round(best_row_bic$GAIC2, 2), "\n")
cat("GD       :", round(best_row_bic$GlobalDeviance, 2), "\n")

cat("\n==================== Best model by AIC ====================\n")
cat("Family   :", best_row_aic$Family, "\n")
cat("df_mu    :", best_row_aic$DF_Mu, "\n")
cat("df_sigma :", best_row_aic$DF_Sigma, "\n")
cat("AIC      :", round(best_row_aic$AIC, 2), "\n")
cat("BIC      :", round(best_row_aic$BIC, 2), "\n")
cat("GAIC2    :", round(best_row_aic$GAIC2, 2), "\n")
cat("GD       :", round(best_row_aic$GlobalDeviance, 2), "\n")

PlotDir <- file.path(OutDir2, "model_comparison_plots")
if (!dir.exists(PlotDir)) dir.create(PlotDir, recursive = TRUE)

BaseSize   <- 14
TitleSize  <- 16
AxisSize   <- 12
LegendSize <- 11
PointSize  <- 3
LineWidth  <- 0.6

theme_model <- function() {
  theme_bw(base_size = BaseSize, base_family = "Arial") +
    theme(
      plot.title   = element_text(size = TitleSize, face = "bold", hjust = 0.5),
      axis.title   = element_text(size = TitleSize - 1),
      axis.text    = element_text(size = AxisSize, color = "black"),
      legend.title = element_text(size = LegendSize),
      legend.text  = element_text(size = LegendSize - 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey85")
    )
}

save_plot <- function(p, filename, w = 10, h = 7, dpi = 600) {
  ggsave(
    filename = file.path(PlotDir, filename),
    plot = p,
    width = w,
    height = h,
    dpi = dpi
  )
}

plot_df <- results_df %>%
  filter(Converged)

if (nrow(plot_df) == 0) {
  warning("No converged models available for plotting.")
} else {
  
  plot_df <- plot_df %>%
    mutate(
      ModelLabel = paste0(Family, " | mu=", DF_Mu, " | sigma=", DF_Sigma)
    )
  
  top_bic <- plot_df %>%
    arrange(BIC, AIC) %>%
    slice(1:20) %>%
    mutate(ModelLabel = factor(ModelLabel, levels = rev(ModelLabel)))
  
  p_bic_top <- ggplot(top_bic, aes(x = ModelLabel, y = BIC, fill = Family)) +
    geom_col(color = "black", linewidth = 0.2) +
    coord_flip() +
    labs(
      title = "Top 20 models ranked by BIC",
      x = NULL,
      y = "BIC"
    ) +
    theme_model()
  
  save_plot(p_bic_top, "01_top20_models_by_BIC.png", w = 11, h = 8)
  
  top_aic <- plot_df %>%
    arrange(AIC, BIC) %>%
    slice(1:20) %>%
    mutate(ModelLabel = factor(ModelLabel, levels = rev(ModelLabel)))
  
  p_aic_top <- ggplot(top_aic, aes(x = ModelLabel, y = AIC, fill = Family)) +
    geom_col(color = "black", linewidth = 0.2) +
    coord_flip() +
    labs(
      title = "Top 20 models ranked by AIC",
      x = NULL,
      y = "AIC"
    ) +
    theme_model()
  
  save_plot(p_aic_top, "02_top20_models_by_AIC.png", w = 11, h = 8)
  
  fam_best_bic <- plot_df %>%
    group_by(Family) %>%
    slice_min(order_by = BIC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Label = paste0("mu=", DF_Mu, ", sigma=", DF_Sigma),
      Family = factor(Family, levels = Family[order(BIC, decreasing = TRUE)])
    )
  
  p_fam_bic <- ggplot(fam_best_bic, aes(x = Family, y = BIC, fill = Family)) +
    geom_col(color = "black", linewidth = 0.2) +
    geom_text(aes(label = Label), hjust = -0.05, size = 3.8, family = "Arial") +
    coord_flip(clip = "off") +
    labs(
      title = "Best BIC model within each family",
      x = NULL,
      y = "BIC"
    ) +
    theme_model() +
    theme(legend.position = "none",
          plot.margin = margin(5.5, 50, 5.5, 5.5))
  
  save_plot(p_fam_bic, "03_best_BIC_within_each_family.png", w = 10, h = 6)
  
  fam_best_aic <- plot_df %>%
    group_by(Family) %>%
    slice_min(order_by = AIC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Label = paste0("mu=", DF_Mu, ", sigma=", DF_Sigma),
      Family = factor(Family, levels = Family[order(AIC, decreasing = TRUE)])
    )
  
  p_fam_aic <- ggplot(fam_best_aic, aes(x = Family, y = AIC, fill = Family)) +
    geom_col(color = "black", linewidth = 0.2) +
    geom_text(aes(label = Label), hjust = -0.05, size = 3.8, family = "Arial") +
    coord_flip(clip = "off") +
    labs(
      title = "Best AIC model within each family",
      x = NULL,
      y = "AIC"
    ) +
    theme_model() +
    theme(legend.position = "none",
          plot.margin = margin(5.5, 50, 5.5, 5.5))
  
  save_plot(p_fam_aic, "04_best_AIC_within_each_family.png", w = 10, h = 6)
  
  for (fam in unique(plot_df$Family)) {
    fam_df <- plot_df %>%
      filter(Family == fam)
    
    p_heat_bic <- ggplot(fam_df, aes(x = factor(DF_Mu), y = factor(DF_Sigma), fill = BIC)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = round(BIC, 1)), size = 3.5, family = "Arial") +
      labs(
        title = paste0("BIC heatmap: ", fam),
        x = expression(df[mu]),
        y = expression(df[sigma]),
        fill = "BIC"
      ) +
      theme_model()
    
    save_plot(p_heat_bic, paste0("05_BIC_heatmap_", fam, ".png"), w = 8, h = 6)
  }
  
  for (fam in unique(plot_df$Family)) {
    fam_df <- plot_df %>%
      filter(Family == fam)
    
    p_heat_aic <- ggplot(fam_df, aes(x = factor(DF_Mu), y = factor(DF_Sigma), fill = AIC)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = round(AIC, 1)), size = 3.5, family = "Arial") +
      labs(
        title = paste0("AIC heatmap: ", fam),
        x = expression(df[mu]),
        y = expression(df[sigma]),
        fill = "AIC"
      ) +
      theme_model()
    
    save_plot(p_heat_aic, paste0("06_AIC_heatmap_", fam, ".png"), w = 8, h = 6)
  }
  
  p_scatter <- ggplot(plot_df, aes(x = AIC, y = BIC, color = Family)) +
    geom_point(size = PointSize, alpha = 0.85) +
    labs(
      title = "AIC vs BIC across all converged models",
      x = "AIC",
      y = "BIC"
    ) +
    theme_model()
  
  save_plot(p_scatter, "07_AIC_vs_BIC_scatter.png", w = 8, h = 6)
  
  family_summary <- plot_df %>%
    group_by(Family) %>%
    summarise(
      N_Converged = n(),
      Best_AIC = min(AIC, na.rm = TRUE),
      Best_BIC = min(BIC, na.rm = TRUE),
      Best_GAIC2 = min(GAIC2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Best_BIC)
  
  write.csv(
    family_summary,
    file.path(OutDir2, "family_summary_AIC_BIC.csv"),
    row.names = FALSE
  )
  
  cat("\nPlots saved to:\n", PlotDir, "\n")
}
