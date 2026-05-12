
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(broom)
  library(stringr)
  library(ggplot2)
  library(ggrain)
})

file1 <- "data/trust_neurodegenerative_phenotypes.xlsx"
outdir  <- "outputs/clinical_dementia"
PlotDir <- "outputs/clinical_dementia/plots"
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

base_needed <- c("z_score", "APOE_Code", "HTN", "DM", "HLP",
                 "MMSE", "MOCA", "AFT", "SDMT", "NPI")
missing_base <- setdiff(base_needed, names(dat))
if(length(missing_base) > 0){
  stop("Missing columns: ", paste(missing_base, collapse = ", "))
}
if(!("BMI_Code" %in% names(dat)) && !("BMI" %in% names(dat))){
  stop("Missing BMI information: provide either BMI_Code or BMI.")
}

dat0 <- dat %>%
  mutate(
    z_score  = as.numeric(z_score),
    APOE_Code = as.numeric(as.character(APOE_Code)),
    MMSE = as.numeric(MMSE),
    MOCA = as.numeric(MOCA),
    AFT  = as.numeric(AFT),
    SDMT = as.numeric(SDMT),
    NPI  = as.numeric(NPI),
    HTN_num = to_num01(HTN),
    DM_num  = to_num01(DM),
    HLP_num = to_num01(HLP),
    BMI_Code = if("BMI_Code" %in% names(.)) {
      to_num01(.data$BMI_Code)
    } else {
      ifelse(as.numeric(BMI) > 25, 1, 0)
    },
    VRS = HTN_num + DM_num + HLP_num + BMI_Code
  ) %>%
  select(z_score, APOE_Code, VRS, MMSE, MOCA, AFT, SDMT, NPI)

run_apoe_linear <- function(data){
  dsub <- data %>%
    select(z_score, APOE_Code) %>%
    filter(complete.cases(.))

  if(nrow(dsub) < 5){
    return(tibble(
      Outcome = "z_score", Predictor = "APOE_Code", Model = "zscore_APOE_linear",
      term = "APOE_Code", beta = NA_real_, std.error = NA_real_, statistic = NA_real_,
      p_value = NA_real_, N = nrow(dsub), R2 = NA_real_, Adj_R2 = NA_real_,
      note = "Too few complete cases"
    ))
  }

  fit <- lm(z_score ~ APOE_Code, data = dsub)
  glance_fit <- glance(fit)
  tidy(fit) %>%
    filter(term == "APOE_Code") %>%
    transmute(
      Outcome = "z_score",
      Predictor = "APOE_Code",
      Model = "zscore_APOE_linear",
      term = term,
      beta = estimate,
      std.error = std.error,
      statistic = statistic,
      p_value = p.value,
      N = nrow(dsub),
      R2 = glance_fit$r.squared,
      Adj_R2 = glance_fit$adj.r.squared,
      note = NA_character_
    )
}

run_cognitive_lm <- function(data, outcome){
  dsub <- data %>%
    select(all_of(c(outcome, "z_score", "VRS"))) %>%
    filter(complete.cases(.))

  if(nrow(dsub) < 5){
    return(tibble(
      Outcome = outcome, Predictor = "z_score", Model = "cognitive_score_zscore_plus_VRS",
      term = "z_score", beta = NA_real_, std.error = NA_real_, statistic = NA_real_,
      p_value = NA_real_, N = nrow(dsub), R2 = NA_real_, Adj_R2 = NA_real_,
      note = "Too few complete cases"
    ))
  }

  fit <- lm(as.formula(paste(outcome, "~ z_score + VRS")), data = dsub)
  glance_fit <- glance(fit)
  tidy(fit) %>%
    filter(term == "z_score") %>%
    transmute(
      Outcome = outcome,
      Predictor = "z_score",
      Model = "cognitive_score_zscore_plus_VRS",
      term = term,
      beta = estimate,
      std.error = std.error,
      statistic = statistic,
      p_value = p.value,
      N = nrow(dsub),
      R2 = glance_fit$r.squared,
      Adj_R2 = glance_fit$adj.r.squared,
      note = NA_character_
    )
}

apoe_result <- run_apoe_linear(dat0)
cognitive_vars <- c("MMSE", "MOCA", "AFT", "SDMT", "NPI")
cognitive_results <- bind_rows(lapply(cognitive_vars, function(v) run_cognitive_lm(dat0, v))) %>%
  mutate(q_value = p.adjust(p_value, method = "fdr"))

write.csv(apoe_result, file.path(outdir, "dementia_minimal_APOE_linear_on_zscore.csv"), row.names = FALSE)
write.csv(cognitive_results, file.path(outdir, "dementia_minimal_cognition_on_zscore_plus_VRS.csv"), row.names = FALSE)
print(apoe_result)
print(cognitive_results)

BaseFamily <- "Arial"
TitleSize <- 30; AxisSize <- 28; TextSize <- 28
LineWidth <- 1.6; AxisWidth <- 1.6; PointSize <- 6
FigW <- 8; FigH <- 6; DPI <- 300; Zthr <- 1.96

theme_my <- function(){
  theme_classic(base_family = BaseFamily) +
    theme(
      plot.title = element_text(size = TitleSize, hjust = 0.5),
      plot.subtitle = element_text(size = TextSize, hjust = 0.5),
      axis.title = element_text(size = AxisSize),
      axis.text = element_text(size = TextSize, colour = "black"),
      axis.line = element_line(linewidth = AxisWidth, colour = "black"),
      axis.ticks = element_line(linewidth = AxisWidth, colour = "black"),
      legend.position = "none",
      panel.grid = element_blank()
    )
}

plot_apoe <- function(data){
  dsub <- data %>% select(z_score, APOE_Code) %>% filter(complete.cases(.))
  dsub$APOE_Code <- factor(dsub$APOE_Code)
  levs <- levels(dsub$APOE_Code)
  cols <- setNames(rep(c("skyblue2", "mediumpurple3", "lightcoral"), length.out = length(levs)), levs)
  p <- ggplot(dsub, aes(x = APOE_Code, y = z_score, fill = APOE_Code)) +
    ggrain::geom_rain(
      jitter_side = "right", point_size = PointSize, alpha = 0.65,
      violin.args = list(alpha = 0.45, linewidth = LineWidth),
      boxplot.args = list(width = 0.15, outlier.shape = NA)
    ) +
    scale_fill_manual(values = cols) +
    geom_hline(yintercept = 0, linewidth = LineWidth, color = "gray55") +
    geom_hline(yintercept = c(-Zthr, Zthr), linetype = "dashed", linewidth = LineWidth, color = "gray55") +
    labs(title = "OEF z-scores by APOE ε4 allele count",
         subtitle = format_p_sub(apoe_result$p_value[1]),
         x = "APOE ε4 allele count", y = "z-score") +
    theme_my()
  ggsave(file.path(PlotDir, "01_zscore_by_APOE_Code.png"), p, width = FigW, height = FigH, dpi = DPI)
  p
}

plot_cognitive <- function(data, v){
  dsub <- data %>% select(all_of(c("z_score", v))) %>% filter(complete.cases(.))
  pval <- cognitive_results$p_value[cognitive_results$Outcome == v][1]
  p <- ggplot(dsub, aes_string(x = "z_score", y = v)) +
    geom_point(size = PointSize, alpha = 0.75, color = "mediumpurple3") +
    geom_smooth(method = "lm", se = TRUE, color = "mediumpurple3", linewidth = LineWidth) +
    labs(title = paste(v, "vs OEF z-scores"),
         subtitle = format_p_sub(pval),
         x = "z-score", y = v) +
    theme_my()
  ggsave(file.path(PlotDir, paste0("02_", v, "_vs_zscore.png")), p, width = FigW, height = FigH, dpi = DPI)
  p
}

plot_apoe(dat0)
invisible(lapply(cognitive_vars, function(v) plot_cognitive(dat0, v)))

cat("\nMinimal dementia outputs saved to:\n", outdir, "\n")
cat("Minimal dementia plots saved to:\n", PlotDir, "\n")
