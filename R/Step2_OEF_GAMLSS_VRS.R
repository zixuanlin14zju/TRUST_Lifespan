
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrain)
})

InFile <- "data/trust_vascular_risk_zscores.csv"
OutDir <- "outputs/vascular_risk"
Phenotype_name <- "OEF"   # used only in titles

Zvar  <- "z_score"
VRS_components <- c("BMI_Code", "HT_Code", "HL_Code", "DB_Code")

RiskVars <- c("BMI_Code", "HT_Code", "HL_Code", "DB_Code", "Smoke", "Drink")
RiskLabels <- c(
  BMI_Code = "BMI Status",
  HT_Code  = "Hypertension",
  HL_Code  = "Hyperlipidemia",
  DB_Code  = "Diabetes",
  Smoke    = "Smoking",
  Drink    = "Alcohol"
)

ExcludeStrokeCol <- "StrokeHistory"  # exclude if == 1; keep 0/NA
ExcludeHeartCol  <- "HeartDisease"   # exclude if == 1; keep 0/NA

FontFamily <- "Arial"
Palette2 <- c("skyblue2", "lightcoral")   # alternating fills
PointColor <- "grey28"
PointAlpha <- 0.45
PointSize  <- 1.1
ViolinAlpha <- 0.55
BoxAlpha    <- 0.85

LineWidth <- 0.6   # ggplot linewidth is in mm-ish; keep modest for print
RefLineColor <- "gray55"

BaseSize  <- 22
AxisSize  <- 24
LabelSize <- 24
TitleSize <- 26

DPI <- 300
FigW <- 7.2   # inches
FigH <- 5.4   # inches

Zrefs <- c(-1.96, 0, 1.96)

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

ensure_outdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

theme_rain <- function() {
  theme_classic(base_family = FontFamily, base_size = BaseSize) +
    theme(
      axis.text  = element_text(size = AxisSize),
      axis.title = element_text(size = LabelSize),
      plot.title = element_text(size = TitleSize)
    )
}

add_z_reference_lines <- function() {
  list(
    geom_hline(yintercept = Zrefs[1], linetype = "dashed", linewidth = LineWidth, color = RefLineColor),
    geom_hline(yintercept = Zrefs[2], linetype = "solid",  linewidth = LineWidth, color = RefLineColor),
    geom_hline(yintercept = Zrefs[3], linetype = "dashed", linewidth = LineWidth, color = RefLineColor)
  )
}

as_factor_clean <- function(x) {
  f <- as.factor(x)
  droplevels(f)
}

save_raincloud_plot <- function(df, x, y, filename,
                                title = NULL, subtitle = NULL,
                                xlab = NULL, ylab = "Z-score",
                                palette = Palette2) {
  
  df[[x]] <- as_factor_clean(df[[x]])
  
  ngrp <- nlevels(df[[x]])
  fills <- rep(palette, length.out = ngrp)
  
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[x]])) +
    geom_rain(
      point.args   = list(size = PointSize, alpha = PointAlpha, color = PointColor),
      boxplot.args = list(width = 0.12, outlier.shape = NA, alpha = BoxAlpha),
      violin.args  = list(alpha = ViolinAlpha, linewidth = LineWidth)
    ) +
    scale_fill_manual(
      name = "Group",
      values = fills
    ) +
    labs(
      title = title %||% "",
      subtitle = subtitle %||% "",
      x = xlab %||% "",
      y = ylab
    ) +
    coord_cartesian(ylim = c(-3, 3)) +
    theme_rain() +
    add_z_reference_lines() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = TitleSize - 2))
  
  ggsave(
    filename = file.path(OutDir, filename),
    plot = p,
    width = FigW, height = FigH, dpi = DPI
  )
}

run_uni_lm <- function(data, outcome, predictor, label = predictor) {
  fml <- as.formula(paste0(outcome, " ~ ", predictor))
  fit <- lm(fml, data = data)
  s <- summary(fit)$coefficients
  rn <- rownames(s)
  idx <- which(rn != "(Intercept)")[1]
  data.frame(
    Variable = label,
    Beta = s[idx, "Estimate"],
    SE = s[idx, "Std. Error"],
    t_value = s[idx, "t value"],
    p_value = s[idx, "Pr(>|t|)"],
    Significant = ifelse(s[idx, "Pr(>|t|)"] < 0.05, "*", ""),
    stringsAsFactors = FALSE
  )
}

format_p <- function(p) {
  if (is.na(p) || !is.finite(p)) return("")
  if (p < 1e-4) {
    return("P < 0.0001")
  } else {
    return(paste0("P = ", format(signif(p, 2), scientific = FALSE, trim = TRUE)))
  }
}

ensure_outdir(OutDir)

dat <- read.csv(InFile)
if (!Zvar %in% names(dat)) stop("Missing z-score column: ", Zvar)

dat$SiteID <- if ("SiteID" %in% names(dat)) as.factor(dat$SiteID) else NULL
dat$Sex    <- if ("Sex"    %in% names(dat)) as.factor(dat$Sex)    else NULL

if ("SiteID" %in% names(dat) && "Sex" %in% names(dat)) {
  cat("\nDescriptive stats of z-score grouped by SiteID and Sex:\n")
  z_summary <- aggregate(
    dat[[Zvar]] ~ SiteID + Sex, data = dat,
    FUN = function(x) {
      x <- x[!is.na(x)]
      c(mean = mean(x), sd = sd(x), median = median(x), n = length(x), se = sd(x)/sqrt(length(x)))
    }
  )
  print(z_summary)
}

dat_filt <- dat
if (ExcludeStrokeCol %in% names(dat_filt)) {
  n0 <- nrow(dat_filt)
  dat_filt <- dat_filt[dat_filt[[ExcludeStrokeCol]] == 0 | is.na(dat_filt[[ExcludeStrokeCol]]), ]
  cat(sprintf("\nExcluded stroke history: %d removed\n", n0 - nrow(dat_filt)))
} else {
  cat("\nStrokeHistory column not found; no exclusion applied.\n")
}

if (ExcludeHeartCol %in% names(dat_filt)) {
  n1 <- nrow(dat_filt)
  dat_filt <- dat_filt[dat_filt[[ExcludeHeartCol]] == 0 | is.na(dat_filt[[ExcludeHeartCol]]), ]
  cat(sprintf("Excluded heart disease history: %d removed\n", n1 - nrow(dat_filt)))
} else {
  cat("HeartDisease column not found; no exclusion applied.\n")
}

cat(sprintf("Final N after exclusions: %d\n", nrow(dat_filt)))

for (v in RiskVars) {
  if (v %in% names(dat_filt)) dat_filt[[v]] <- as.factor(dat_filt[[v]])
}

missing_vrs_cols <- setdiff(VRS_components, names(dat_filt))
if (length(missing_vrs_cols) > 0) {
  stop("Missing VRS component columns: ", paste(missing_vrs_cols, collapse = ", "))
}

for (v in VRS_components) {
  dat_filt[[paste0(v, "_num")]] <- suppressWarnings(as.integer(as.character(dat_filt[[v]])))
  dat_filt[[paste0(v, "_num")]][!(dat_filt[[paste0(v, "_num")]] %in% c(0, 1))] <- NA_integer_
}

dat_filt$VRS_complete <- complete.cases(dat_filt[, paste0(VRS_components, "_num")])

dat_filt$VRS_any_info <- rowSums(!is.na(dat_filt[, paste0(VRS_components, "_num"), drop = FALSE])) > 0

dat_filt$VRS_numeric <- NA_integer_
dat_filt$VRS_numeric[dat_filt$VRS_complete] <- rowSums(
  dat_filt[dat_filt$VRS_complete, paste0(VRS_components, "_num"), drop = FALSE]
)

cat("\nVRS calculated from: BMI_Code + HT_Code + HL_Code + DB_Code\n")
cat(sprintf("Subjects with complete data for all 4 VRS components: %d / %d\n",
            sum(dat_filt$VRS_complete), nrow(dat_filt)))
cat(sprintf("Subjects with at least one available VRS component: %d / %d\n",
            sum(dat_filt$VRS_any_info), nrow(dat_filt)))
cat(sprintf("Subjects with no VRS information at all: %d / %d\n",
            sum(!dat_filt$VRS_any_info), nrow(dat_filt)))

if (all(c("Smoke", "Drink") %in% names(dat_filt))) {
  
  smoke_num <- suppressWarnings(as.integer(as.character(dat_filt$Smoke)))
  drink_num <- suppressWarnings(as.integer(as.character(dat_filt$Drink)))
  
  smoke_num[!(smoke_num %in% c(0, 1))] <- NA_integer_
  drink_num[!(drink_num %in% c(0, 1))] <- NA_integer_
  
  dat_filt$Lifestyle <- smoke_num + drink_num  # 0,1,2 (or NA)
  dat_filt$Lifestyle <- factor(dat_filt$Lifestyle, levels = 0:2, labels = c("0", "1", "2"), ordered = TRUE)
  
} else {
  cat("\nSmoke/Drink not found; Lifestyle score will not be created.\n")
}

cat("\nDescriptive stats after exclusions:\n")
desc_stats <- data.frame(
  Variable = character(),
  N = integer(),
  Summary = character(),
  stringsAsFactors = FALSE
)

desc_stats <- rbind(desc_stats,
                    data.frame(
                      Variable = "Total sample after exclusions",
                      N = nrow(dat_filt),
                      Summary = "100%"
                    ),
                    data.frame(
                      Variable = "Complete data for all 4 VRS components",
                      N = sum(dat_filt$VRS_complete),
                      Summary = sprintf("%.1f%%", 100 * mean(dat_filt$VRS_complete))
                    ),
                    data.frame(
                      Variable = "Incomplete data for >=1 VRS component",
                      N = sum(dat_filt$VRS_any_info),
                      Summary = sprintf("%.1f%%", 100 * mean(!dat_filt$VRS_complete))
                    )
)

for (v in VRS_components) {
  v_num <- dat_filt[[paste0(v, "_num")]]
  desc_stats <- rbind(desc_stats,
                      data.frame(
                        Variable = paste0(v, " available"),
                        N = sum(!is.na(v_num)),
                        Summary = sprintf("%.1f%%", 100 * mean(!is.na(v_num)))
                      )
  )
}

if ("VRS_numeric" %in% names(dat_filt)) {
  vrs_nonmiss <- dat_filt$VRS_numeric[is.finite(dat_filt$VRS_numeric)]
  
  desc_stats <- rbind(desc_stats,
                      data.frame(
                        Variable = "VRS (complete cases only)",
                        N = length(vrs_nonmiss),
                        Summary = sprintf("%.2f \u00B1 %.2f",
                                          mean(vrs_nonmiss, na.rm = TRUE),
                                          sd(vrs_nonmiss, na.rm = TRUE))
                      )
  )
  
  vrs_tab <- table(factor(vrs_nonmiss, levels = 0:4), useNA = "no")
  vrs_denom <- sum(vrs_tab)
  
  for (lvl in names(vrs_tab)) {
    desc_stats <- rbind(desc_stats,
                        data.frame(
                          Variable = paste0("VRS (", lvl, ")"),
                          N = as.integer(vrs_tab[[lvl]]),
                          Summary = sprintf("%.1f%%", 100 * vrs_tab[[lvl]] / vrs_denom)
                        )
    )
  }
}

if ("BMI_Code_num" %in% names(dat_filt)) {
  bmi_num <- dat_filt$BMI_Code_num
  desc_stats <- rbind(desc_stats, data.frame(
    Variable = "BMI_Code",
    N = sum(is.finite(bmi_num)),
    Summary = sprintf("%.2f \u00B1 %.2f",
                      mean(bmi_num, na.rm = TRUE),
                      sd(bmi_num, na.rm = TRUE))
  ))
}

for (v in RiskVars) {
  if (v %in% names(dat_filt)) {
    tab <- table(dat_filt[[v]], useNA = "no")
    denom <- sum(tab)
    for (lvl in names(tab)) {
      if (!is.na(lvl)) {
        desc_stats <- rbind(desc_stats, data.frame(
          Variable = paste0(v, " (", lvl, ")"),
          N = as.integer(tab[[lvl]]),
          Summary = sprintf("%.1f%%", 100 * tab[[lvl]] / denom)
        ))
      }
    }
  }
}

if ("Lifestyle" %in% names(dat_filt)) {
  tab <- table(dat_filt$Lifestyle, useNA = "ifany")
  denom <- sum(tab)
  
  for (lvl in names(tab)) {
    if (!is.na(lvl)) {
      desc_stats <- rbind(desc_stats, data.frame(
        Variable = paste0("Lifestyle (", lvl, ")"),
        N = as.integer(tab[[lvl]]),
        Summary = sprintf("%.1f%%", 100 * tab[[lvl]] / denom),
        stringsAsFactors = FALSE
      ))
    }
  }
}

print(desc_stats)

cat("\nUnivariate linear models (one predictor at a time):\n")

results_df <- data.frame(
  Variable = character(),
  Beta = numeric(),
  SE = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

if ("VRS_numeric" %in% names(dat_filt)) {
  dat_vrs <- dat_filt[is.finite(dat_filt$VRS_numeric) & is.finite(dat_filt[[Zvar]]), ]
  if (nrow(dat_vrs) >= 5) {
    results_df <- rbind(results_df, run_uni_lm(dat_vrs, Zvar, "VRS_numeric", "VRS"))
  }
}

for (v in RiskVars) {
  if (v %in% names(dat_filt)) {
    results_df <- rbind(results_df, run_uni_lm(dat_filt, Zvar, v, v))
  }
}

if ("Lifestyle" %in% names(dat_filt)) {
  dat_filt$Lifestyle_numeric <- suppressWarnings(as.numeric(as.character(dat_filt$Lifestyle)))
  results_df <- rbind(results_df, run_uni_lm(dat_filt, Zvar, "Lifestyle_numeric", "Lifestyle (0–2)"))
}

print(results_df)

cat("\nMultivariable linear model with all VRS components entered simultaneously:\n")

vrs_component_num_vars <- paste0(VRS_components, "_num")

dat_multi_vrs <- dat_filt[
  complete.cases(dat_filt[, c(Zvar, vrs_component_num_vars)]) &
    is.finite(dat_filt[[Zvar]]),
]

cat(sprintf("N for multivariable VRS-component model: %d\n", nrow(dat_multi_vrs)))

multivariable_vrs_results <- data.frame(
  Variable = character(),
  Beta = numeric(),
  SE = numeric(),
  t_value = numeric(),
  p_value = numeric(),
  Significant = character(),
  stringsAsFactors = FALSE
)

if (nrow(dat_multi_vrs) >= 10) {
  
  fml_multi_vrs <- as.formula(
    paste0(Zvar, " ~ ", paste(vrs_component_num_vars, collapse = " + "))
  )
  
  fit_multi_vrs <- lm(fml_multi_vrs, data = dat_multi_vrs)
  print(summary(fit_multi_vrs))
  
  coef_tab <- summary(fit_multi_vrs)$coefficients
  coef_tab <- coef_tab[rownames(coef_tab) != "(Intercept)", , drop = FALSE]
  
  multivariable_vrs_results <- data.frame(
    Variable = rownames(coef_tab),
    Beta = coef_tab[, "Estimate"],
    SE = coef_tab[, "Std. Error"],
    t_value = coef_tab[, "t value"],
    p_value = coef_tab[, "Pr(>|t|)"],
    Significant = ifelse(coef_tab[, "Pr(>|t|)"] < 0.05, "*", ""),
    stringsAsFactors = FALSE
  )
  
  multivariable_vrs_results$Label <- dplyr::recode(
    multivariable_vrs_results$Variable,
    BMI_Code_num = "BMI Status",
    HT_Code_num  = "Hypertension",
    HL_Code_num  = "Hyperlipidemia",
    DB_Code_num  = "Diabetes"
  )
  
  dat_multi_scaled <- dat_multi_vrs
  dat_multi_scaled[[Zvar]] <- scale(dat_multi_scaled[[Zvar]])[, 1]
  for (v in vrs_component_num_vars) {
    dat_multi_scaled[[v]] <- scale(dat_multi_scaled[[v]])[, 1]
  }
  
  fit_multi_vrs_std <- lm(fml_multi_vrs, data = dat_multi_scaled)
  coef_tab_std <- summary(fit_multi_vrs_std)$coefficients
  coef_tab_std <- coef_tab_std[rownames(coef_tab_std) != "(Intercept)", , drop = FALSE]
  
  multivariable_vrs_results$Std_Beta <- coef_tab_std[
    match(multivariable_vrs_results$Variable, rownames(coef_tab_std)),
    "Estimate"
  ]
  
  multivariable_vrs_results <- multivariable_vrs_results[
    order(abs(multivariable_vrs_results$Std_Beta), decreasing = TRUE),
  ]
  
  cat("\nMultivariable VRS-component results ordered by absolute standardized beta:\n")
  print(multivariable_vrs_results)
  
  write.csv(
    multivariable_vrs_results,
    file = file.path(OutDir, "multivariable_VRS_component_model_results.csv"),
    row.names = FALSE
  )
  
} else {
  cat("Not enough complete cases for multivariable VRS-component model.\n")
}

if (exists("multivariable_vrs_results") &&
    nrow(multivariable_vrs_results) > 0 &&
    "Std_Beta" %in% names(multivariable_vrs_results)) {
  
  risk_colors <- c(
    "BMI Status" = "skyblue2",
    "Hypertension" = "lightcoral",
    "Hyperlipidemia" = "mediumpurple3",
    "Diabetes" = "goldenrod3"
  )
  
  plot_df <- multivariable_vrs_results
  
  plot_df$Label <- factor(
    plot_df$Label,
    levels = plot_df$Label[order(plot_df$Std_Beta)]
  )
  
  p_beta <- ggplot(
    plot_df,
    aes(
      x = Label,
      y = Std_Beta,
      fill = Label
    )
  ) +
    geom_col(width = 0.68) +
    geom_hline(yintercept = 0, linewidth = LineWidth, color = RefLineColor) +
    coord_flip() +
    scale_fill_manual(values = risk_colors) +
    labs(
      title = "Relative contribution of vascular risk components",
      x = "",
      y = "Standardized beta"
    ) +
    theme_classic(base_family = FontFamily, base_size = BaseSize) +
    theme(
      axis.text  = element_text(size = AxisSize),
      axis.title = element_text(size = LabelSize),
      plot.title = element_text(size = TitleSize - 4, hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(t = 12, r = 35, b = 12, l = 12)
    )
  
  ggsave(
    filename = file.path(OutDir, "multivariable_VRS_component_standardized_beta.png"),
    plot = p_beta,
    width = 9,      # wider than before
    height = 5.4,
    dpi = DPI,
    limitsize = FALSE
  )
}

cat("\nSensitivity analysis: VRS treated as categorical variable:\n")

anova_vrs_results <- data.frame(
  Model = character(),
  Df = numeric(),
  F_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

if ("VRS_numeric" %in% names(dat_filt)) {
  
  dat_vrs_cat <- dat_filt[
    is.finite(dat_filt$VRS_numeric) & is.finite(dat_filt[[Zvar]]),
  ]
  
  dat_vrs_cat$VRS_cat <- factor(
    as.integer(round(dat_vrs_cat$VRS_numeric)),
    levels = 0:4
  )
  
  dat_vrs_cat$VRS_cat <- droplevels(dat_vrs_cat$VRS_cat)
  
  if (nrow(dat_vrs_cat) >= 5 && nlevels(dat_vrs_cat$VRS_cat) >= 2) {
    
    fit_vrs_cat <- lm(as.formula(paste0(Zvar, " ~ VRS_cat")), data = dat_vrs_cat)
    aov_vrs_cat <- anova(fit_vrs_cat)
    
    print(aov_vrs_cat)
    
    anova_vrs_results <- rbind(
      anova_vrs_results,
      data.frame(
        Model = "VRS categorical 0-4",
        Df = aov_vrs_cat["VRS_cat", "Df"],
        F_value = aov_vrs_cat["VRS_cat", "F value"],
        p_value = aov_vrs_cat["VRS_cat", "Pr(>F)"],
        stringsAsFactors = FALSE
      )
    )
    
    cat("\nGroup-wise z-score summary by VRS category:\n")
    vrs_cat_summary <- aggregate(
      dat_vrs_cat[[Zvar]] ~ VRS_cat,
      data = dat_vrs_cat,
      FUN = function(x) {
        x <- x[is.finite(x)]
        c(
          n = length(x),
          mean = mean(x),
          sd = sd(x),
          se = sd(x) / sqrt(length(x))
        )
      }
    )
    print(vrs_cat_summary)
  }
}

print(anova_vrs_results)

for (v in RiskVars) {
  
  if (!v %in% names(dat_filt)) next
  
  label <- RiskLabels[v]
  
  dfp <- dat_filt[, c(Zvar, v)]
  names(dfp) <- c("z", "grp")
  dfp <- dfp[is.finite(dfp$z) & !is.na(dfp$grp), ]
  
  if (nrow(dfp) < 5 || nlevels(as.factor(dfp$grp)) < 2) next
  
  p_sub <- ""
  if (v %in% results_df$Variable) {
    p_sub <- format_p(results_df$p_value[match(v, results_df$Variable)])
  }
  
  save_raincloud_plot(
    df = transform(dfp, grp = as.factor(grp)),
    x  = "grp",
    y  = "z",
    filename = paste0("zscore_", v, "_raincloud.png"),
    title = paste0(Phenotype_name, " z-scores by ", label),
    subtitle = p_sub,
    xlab = label,
    ylab = "z-score"
  )
}

p_sub_vrs <- ""
if ("VRS" %in% results_df$Variable) {
  p_sub_vrs <- format_p(results_df$p_value[match("VRS", results_df$Variable)])
}

if ("VRS_numeric" %in% names(dat_filt)) {
  
  dfv <- dat_filt[, c(Zvar, "VRS_numeric", "VRS_complete")]
  names(dfv) <- c("z", "vrs", "complete")
  dfv <- dfv[dfv$complete & is.finite(dfv$z) & is.finite(dfv$vrs), ]
  
  dfv$vrs_int <- as.integer(round(dfv$vrs))
  dfv <- dfv[dfv$vrs_int %in% 0:4, ]
  
  dfv$VRS_group <- factor(dfv$vrs_int, levels = 0:4, labels = as.character(0:4), ordered = TRUE)
  
  if (nrow(dfv) >= 5 && nlevels(dfv$VRS_group) >= 2) {
    
    vrs_purples <- colorRampPalette(c("thistle3", "mediumpurple3", "purple4"))(5)
    names(vrs_purples) <- as.character(0:4)
    
    p <- ggplot(dfv, aes(x = VRS_group, y = z, fill = VRS_group)) +
      geom_rain(
        point.args   = list(size = PointSize, alpha = PointAlpha, color = PointColor),
        boxplot.args = list(width = 0.12, outlier.shape = NA, alpha = BoxAlpha),
        violin.args  = list(alpha = ViolinAlpha, linewidth = LineWidth)
      ) +
      scale_fill_manual(
        name = "VRS",
        values = vrs_purples
      ) +
      scale_x_discrete(drop = FALSE) +
      labs(
        title = paste0(Phenotype_name, " z-scores across VRS"),
        subtitle = p_sub_vrs,
        x = "VRS",
        y = "z-score"
      ) +
      coord_cartesian(ylim = c(-3, 3)) +
      theme_rain() +
      add_z_reference_lines() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = TitleSize - 2))
    
    ggsave(
      filename = file.path(OutDir, "zscore_vs_VRS_raincloud.png"),
      plot = p,
      width = FigW, height = FigH, dpi = DPI
    )
  }
}

p_sub_ls <- ""
if ("Lifestyle (0–2)" %in% results_df$Variable) {
  p_sub_ls <- format_p(results_df$p_value[match("Lifestyle (0–2)", results_df$Variable)])
}

if ("Lifestyle" %in% names(dat_filt)) {
  
  dfl <- dat_filt[, c(Zvar, "Lifestyle")]
  names(dfl) <- c("z", "ls")
  dfl <- dfl[is.finite(dfl$z) & !is.na(dfl$ls), ]
  
  if (nrow(dfl) >= 5 && nlevels(as.factor(dfl$ls)) >= 2) {
    
    ls_blues <- colorRampPalette(c("palegreen3", "mediumseagreen", "seagreen"))(3)
    names(ls_blues) <- c("0", "1", "2")
    
    p_ls <- ggplot(dfl, aes(x = ls, y = z, fill = ls)) +
      geom_rain(
        point.args   = list(size = PointSize, alpha = PointAlpha, color = PointColor),
        boxplot.args = list(width = 0.12, outlier.shape = NA, alpha = BoxAlpha),
        violin.args  = list(alpha = ViolinAlpha, linewidth = LineWidth)
      ) +
      scale_fill_manual(
        name = "Lifestyle",
        values = ls_blues
      ) +
      scale_x_discrete(drop = FALSE) +
      labs(
        title = paste0(Phenotype_name, " z-scores across Lifestyle"),
        subtitle = p_sub_ls,
        x = "Lifestyle (Smoking + Alcohol)",
        y = "z-score"
      ) +
      coord_cartesian(ylim = c(-3, 3)) +
      theme_rain() +
      add_z_reference_lines() +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = TitleSize - 4))
    
    ggsave(
      filename = file.path(OutDir, "zscore_vs_Lifestyle_raincloud.png"),
      plot = p_ls,
      width = FigW, height = FigH, dpi = DPI
    )
  }
}
