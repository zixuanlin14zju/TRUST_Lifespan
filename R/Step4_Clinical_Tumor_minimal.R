
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(broom)
  library(stringr)
  library(ggplot2)
  library(ggrain)
})

file1 <- "data/trust_tumor_phenotypes.xlsx"
outdir  <- "outputs/clinical_tumor"
PlotDir <- "outputs/clinical_tumor/plots"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if(!dir.exists(PlotDir)) dir.create(PlotDir, recursive = TRUE)

dat <- read_excel(file1, sheet = "Sheet1")

clean_names_simple <- function(x){
  x <- gsub("\\s+", "_", x)
  x <- gsub("-", "_", x)
  x <- gsub("/", "_", x)
  x <- gsub("\\(", "", x)
  x <- gsub("\\)", "", x)
  x <- gsub("\\.", "_", x)
  x <- gsub("β", "beta", x)
  x <- gsub("τ", "tau", x)
  x <- gsub("α", "alpha", x)
  x <- gsub("Aβ", "Abeta", x)
  x <- gsub("[^[:alnum:]_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

names(dat) <- clean_names_simple(names(dat))

to_num01 <- function(x){
  suppressWarnings(as.numeric(as.character(x)))
}

format_p_sub <- function(p){
  if(is.na(p) || !is.finite(p)) return("")
  if(p < 1e-4) "P < 0.0001" else paste0("P = ", format(signif(p, 2), scientific = FALSE, trim = TRUE))
}

base_needed <- c("z_score", "Age", "Sex", "HTN", "DM", "HLP",
                 "Tumor_Grade", "Ki67con", "Ki67dis", "IDH1")
missing_base <- setdiff(base_needed, names(dat))
if(length(missing_base) > 0){
  stop("Missing columns: ", paste(missing_base, collapse = ", "))
}
if(!("BMI_Code" %in% names(dat)) && !("BMI" %in% names(dat))){
  stop("Missing BMI information: provide either BMI_Code or BMI.")
}

dat0 <- dat %>%
  mutate(
    z_score     = as.numeric(z_score),
    Age         = as.numeric(Age),
    Sex         = factor(Sex),
    Tumor_Grade = as.numeric(Tumor_Grade),
    Ki67con     = as.numeric(Ki67con),
    Ki67dis     = as.numeric(Ki67dis),
    IDH1        = factor(IDH1),
    HTN_num     = to_num01(HTN),
    DM_num      = to_num01(DM),
    HLP_num     = to_num01(HLP),
    BMI_Code    = if("BMI_Code" %in% names(.)) {
      to_num01(.data$BMI_Code)
    } else {
      ifelse(as.numeric(BMI) > 25, 1, 0)
    },
    VRS         = HTN_num + DM_num + HLP_num + BMI_Code
  ) %>%
  select(z_score, Age, Sex, VRS, Tumor_Grade, Ki67con, Ki67dis, IDH1)

dat0 <- dat0 %>%
  mutate(
    Ki67dis = ifelse(is.na(Ki67dis) & !is.na(Ki67con) & abs(Ki67con - 0.2) < 1e-8, 1, Ki67dis)
  )

run_lm_global <- function(data, predictor, covariates = c("Age", "Sex", "VRS")){
  vars_use <- unique(c("z_score", predictor, covariates))
  dsub <- data %>%
    select(all_of(vars_use)) %>%
    filter(complete.cases(.)) %>%
    droplevels()

  if(nrow(dsub) < 5){
    return(tibble(
      Predictor = predictor, term = predictor, beta = NA_real_, std.error = NA_real_,
      statistic = NA_real_, p_coef = NA_real_, p_global = NA_real_,
      N = nrow(dsub), R2 = NA_real_, Adj_R2 = NA_real_, note = "Too few complete cases"
    ))
  }

  if(is.factor(dsub[[predictor]]) && nlevels(dsub[[predictor]]) < 2){
    return(tibble(
      Predictor = predictor, term = predictor, beta = NA_real_, std.error = NA_real_,
      statistic = NA_real_, p_coef = NA_real_, p_global = NA_real_,
      N = nrow(dsub), R2 = NA_real_, Adj_R2 = NA_real_, note = "Predictor has <2 levels"
    ))
  }

  fit <- lm(as.formula(paste("z_score ~", paste(c(predictor, covariates), collapse = " + "))), data = dsub)
  fit_glance <- glance(fit)
  fit_tidy <- tidy(fit)
  fit_pred <- fit_tidy %>%
    filter(term == predictor | str_detect(term, paste0("^", predictor))) %>%
    rename(beta = estimate, p_coef = p.value)

  drop_res <- tryCatch(drop1(fit, test = "F"), error = function(e) NULL)
  p_global <- NA_real_
  if(!is.null(drop_res) && predictor %in% rownames(drop_res)){
    p_global <- drop_res[predictor, "Pr(>F)"]
  }

  fit_pred %>%
    mutate(
      Predictor = predictor,
      Model = "zscore_predictor_plus_Age_Sex_VRS",
      N = nrow(dsub),
      R2 = fit_glance$r.squared,
      Adj_R2 = fit_glance$adj.r.squared,
      p_global = p_global,
      note = NA_character_
    ) %>%
    select(Predictor, Model, term, beta, std.error, statistic, p_coef, p_global, N, R2, Adj_R2, note)
}

predictors <- c("Tumor_Grade", "Ki67con", "Ki67dis", "IDH1")
results <- bind_rows(lapply(predictors, function(v) run_lm_global(dat0, v))) %>%
  mutate(q_global = p.adjust(p_global, method = "fdr"),
         q_coef   = p.adjust(p_coef, method = "fdr"))

summary_table <- results %>%
  group_by(Predictor, Model) %>%
  summarise(
    N = unique(N),
    R2 = unique(R2),
    Adj_R2 = unique(Adj_R2),
    p_global = unique(p_global),
    p_coef_min = suppressWarnings(min(p_coef, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    p_coef_min = ifelse(is.finite(p_coef_min), p_coef_min, NA_real_),
    q_global = p.adjust(p_global, method = "fdr"),
    q_coef   = p.adjust(p_coef_min, method = "fdr")
  )

write.csv(results, file.path(outdir, "tumor_minimal_model_terms.csv"), row.names = FALSE)
write.csv(summary_table, file.path(outdir, "tumor_minimal_summary.csv"), row.names = FALSE)
print(summary_table)

get_p <- function(pred){
  p <- summary_table$p_global[summary_table$Predictor == pred][1]
  format_p_sub(p)
}

BaseFamily <- "Arial"
TitleSize <- 32; AxisSize <- 32; TextSize <- 32
LineWidth <- 1.2; AxisWidth <- 1.2; PointSize <- 5
FigW <- 8; FigH <- 6; DPI <- 300; Zthr <- 1.96

theme_my <- function(){
  theme_classic(base_family = BaseFamily) +
    theme(
      plot.title = element_text(size = TitleSize, hjust = 0.5),
      plot.subtitle = element_text(size = TitleSize - 4, hjust = 0.5),
      axis.title = element_text(size = AxisSize),
      axis.text = element_text(size = TextSize, colour = "black"),
      axis.line = element_line(linewidth = AxisWidth, colour = "black"),
      axis.ticks = element_line(linewidth = AxisWidth, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank()
    )
}

plot_raincloud <- function(data, xvar, outname, xlab = xvar, title = NULL, subtitle = NULL){
  dsub <- data %>% select(all_of(c(xvar, "z_score"))) %>% filter(complete.cases(.))
  dsub[[xvar]] <- factor(dsub[[xvar]])
  levs <- levels(dsub[[xvar]])
  cols <- setNames(rep(c("skyblue2", "mediumpurple3", "lightcoral", "tan3", "darkseagreen3"),
                       length.out = length(levs)), levs)

  p <- ggplot(dsub, aes_string(x = xvar, y = "z_score", fill = xvar)) +
    ggrain::geom_rain(
      jitter_side = "right", point_size = PointSize, alpha = 0.65,
      violin.args = list(alpha = 0.45, linewidth = LineWidth),
      boxplot.args = list(width = 0.15, outlier.shape = NA)
    ) +
    scale_fill_manual(values = cols) +
    geom_hline(yintercept = 0, linewidth = LineWidth, color = "gray55") +
    geom_hline(yintercept = c(-Zthr, Zthr), linetype = "dashed", linewidth = LineWidth, color = "gray55") +
    labs(title = title, subtitle = subtitle, x = xlab, y = "z-score") +
    theme_my()
  ggsave(file.path(PlotDir, outname), p, width = FigW, height = FigH, dpi = DPI)
  p
}

plot_scatter <- function(data, xvar, outname, xlab = xvar, title = NULL, subtitle = NULL){
  dsub <- data %>% select(all_of(c(xvar, "z_score"))) %>% filter(complete.cases(.))
  p <- ggplot(dsub, aes_string(x = xvar, y = "z_score")) +
    geom_point(size = PointSize, alpha = 0.75, color = "mediumpurple3") +
    geom_smooth(method = "lm", se = TRUE, color = "mediumpurple3", linewidth = LineWidth) +
    labs(title = title, subtitle = subtitle, x = xlab, y = "z-score") +
    theme_my()
  ggsave(file.path(PlotDir, outname), p, width = FigW, height = FigH, dpi = DPI)
  p
}

plot_dat <- dat0 %>%
  mutate(
    Tumor_Grade = factor(Tumor_Grade),
    Ki67dis = factor(Ki67dis, levels = sort(unique(Ki67dis)), ordered = TRUE),
    IDH1 = factor(IDH1)
  )

if(all(na.omit(unique(as.character(plot_dat$IDH1))) %in% c("0", "1"))){
  plot_dat$IDH1 <- factor(plot_dat$IDH1, levels = c("0", "1"), labels = c("Wildtype", "Mutant"))
}
if(all(na.omit(unique(as.character(plot_dat$Ki67dis))) %in% c("0", "1", "2"))){
  plot_dat$Ki67dis <- factor(plot_dat$Ki67dis, levels = c("0", "1", "2"),
                             labels = c("Low", "Mid", "High"), ordered = TRUE)
}

plot_raincloud(plot_dat, "Tumor_Grade", "01_zscore_by_Tumor_Grade.png",
               xlab = "Tumor grade", title = "OEF z-scores by Tumor Grade", subtitle = get_p("Tumor_Grade"))
plot_raincloud(plot_dat, "Ki67dis", "02_zscore_by_Ki67_Group.png",
               xlab = "Ki-67 group", title = "OEF z-scores by Ki-67 Group", subtitle = get_p("Ki67dis"))
plot_raincloud(plot_dat, "IDH1", "03_zscore_by_IDH1.png",
               xlab = "IDH1", title = "OEF z-scores by IDH1 Status", subtitle = get_p("IDH1"))
plot_scatter(plot_dat, "Ki67con", "04_zscore_vs_Ki67_continuous.png",
             xlab = "Ki-67 (continuous)", title = "OEF z-scores vs Ki-67", subtitle = get_p("Ki67con"))

cat("\nMinimal tumor outputs saved to:\n", outdir, "\n")
cat("Minimal tumor plots saved to:\n", PlotDir, "\n")
