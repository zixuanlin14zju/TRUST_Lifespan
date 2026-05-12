
suppressPackageStartupMessages({
  library(gamlss)
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(ggrain)
})

BestModelPath <- "outputs/lifespan_model/data/best_model.rds"
HCPath        <- "data/trust_lifespan_hc_zscores.csv"
PatientPath   <- "data/trust_disease_patients.csv"

OutDir        <- "outputs/disease_pattern"
if(!dir.exists(OutDir)) dir.create(OutDir, recursive = TRUE)

PhenotypeName <- "OEF"             # for labels only
PhenotypeCol  <- "phenotype"       # patient_data column name that stores OEF values
DiagnosisCol  <- "Diagnosis"
MinSamplesPerDx <- 5

BaseFamily    <- "Arial"
AxisCex       <- 22
TitleCex      <- 24
LineWidth     <- 0.8
PointSize     <- 1.6
FigW          <- 12
FigH          <- 7
DPI           <- 300

Zthr <- 1.96

best_model   <- readRDS(BestModelPath)
M_HC         <- read.csv(HCPath)
patient_data <- read.csv(PatientPath)

stopifnot(DiagnosisCol %in% names(patient_data))
stopifnot(PhenotypeCol %in% names(patient_data))
stopifnot("Age" %in% names(patient_data))
stopifnot("Sex" %in% names(patient_data))

cat("\n=== Patient data overview ===\n")
cat(sprintf("Age range: %.1f - %.1f years\n", min(patient_data$Age, na.rm=TRUE), max(patient_data$Age, na.rm=TRUE)))
cat(sprintf("Sex counts (Sex==1 Female, Sex==0 Male): F=%d, M=%d\n",
            sum(patient_data$Sex == 1, na.rm=TRUE),
            sum(patient_data$Sex == 0, na.rm=TRUE)))
cat("Diagnosis counts:\n")
print(table(patient_data[[DiagnosisCol]]))

if ("SiteID" %in% names(patient_data) && "SiteID" %in% names(M_HC)) {
  train_site_levels <- levels(factor(M_HC$SiteID))
  patient_data$SiteID <- factor(patient_data$SiteID, levels = train_site_levels)
}

need_site <- "SiteID" %in% names(patient_data) || any(grepl("SiteID", deparse(formula(best_model, "mu"))))

newdata_pat <- data.frame(
  Age = patient_data$Age,
  Sex = patient_data$Sex
)

if ("SiteID" %in% names(patient_data) && "SiteID" %in% names(M_HC)) {
  train_site_levels <- levels(factor(M_HC$SiteID))
  patient_data$SiteID <- factor(patient_data$SiteID, levels = train_site_levels)
  
  site_dummy <- train_site_levels[which(!is.na(train_site_levels))[1]]
  newdata_pat$SiteID <- factor(
    rep(site_dummy, nrow(newdata_pat)),
    levels = train_site_levels
  )
}

cat("\n=== Refitting clean GAMLSS model on HC ===\n")
con <- gamlss.control(n.cyc = 200)
df_mu    <- 2
df_sigma <- 1
clean_model <- gamlss(
  phenotype ~ bs(Age, df = df_mu) * Sex + random(as.factor(SiteID)),
  sigma.fo  = ~ bs(Age, df = df_sigma) + Sex,
  nu.fo     = ~ 1,
  tau.fo    = ~ 1,
  family    = BCTo,
  data      = M_HC,
  control   = con
)

cat("\n=== Predicting distribution parameters (population-level; random=zero) ===\n")
pred_all <- predictAll(clean_model, newdata = newdata_pat, random = "zero")

mu_patient    <- pred_all$mu
sigma_patient <- pred_all$sigma
nu_patient    <- pred_all$nu
tau_patient   <- pred_all$tau

cat(sprintf("mu=%d, sigma=%d, nu=%d, tau=%d | Patient rows=%d\n",
            length(mu_patient), length(sigma_patient), length(nu_patient), length(tau_patient),
            nrow(patient_data)))

y <- patient_data[[PhenotypeCol]]
cumulative_prob <- pBCTo(y, mu = mu_patient, sigma = sigma_patient, nu = nu_patient, tau = tau_patient)

cumulative_prob <- pmin(pmax(cumulative_prob, 1e-10), 1 - 1e-10)

patient_data$z_score <- qnorm(cumulative_prob)

cat("\n=== Patient Z-score summary ===\n")
cat(sprintf("Mean: %.3f | SD: %.3f | Range: [%.3f, %.3f]\n",
            mean(patient_data$z_score, na.rm=TRUE),
            sd(patient_data$z_score, na.rm=TRUE),
            min(patient_data$z_score, na.rm=TRUE),
            max(patient_data$z_score, na.rm=TRUE)))

out_z <- file.path(OutDir, "disease_zscores.csv")
write.csv(patient_data, out_z, row.names = FALSE)
cat(sprintf("\nSaved patient Z-scores to: %s\n", out_z))

dx_counts <- table(patient_data[[DiagnosisCol]])
valid_dx  <- names(dx_counts[dx_counts >= MinSamplesPerDx])

cat("\n=== Diagnosis filtering ===\n")
cat(sprintf("Total diagnosis categories: %d\n", length(dx_counts)))
cat(sprintf("Keeping dx with N >= %d: %d categories\n", MinSamplesPerDx, length(valid_dx)))
print(dx_counts)

patient_f <- patient_data %>%
  filter(.data[[DiagnosisCol]] %in% valid_dx) %>%
  mutate(Diagnosis = factor(.data[[DiagnosisCol]], levels = valid_dx))

if (length(valid_dx) < 2) stop("Need at least 2 diagnosis groups (after filtering) for group comparisons.")

desc_stats <- patient_f %>%
  group_by(Diagnosis) %>%
  summarise(
    N = n(),
    Mean_Z = mean(z_score, na.rm = TRUE),
    SD_Z   = sd(z_score, na.rm = TRUE),
    SE_Z   = SD_Z / sqrt(N),
    Median_Z = median(z_score, na.rm = TRUE),
    Min_Z  = min(z_score, na.rm = TRUE),
    Max_Z  = max(z_score, na.rm = TRUE),
    CI95_L = Mean_Z - Zthr * SE_Z,
    CI95_U = Mean_Z + Zthr * SE_Z,
    Abnormal_percent = mean(abs(z_score) > Zthr, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_Z))

cat("\n=== Descriptive stats by diagnosis (filtered) ===\n")
print(desc_stats)

results_list <- list()

for (dx in levels(patient_f$Diagnosis)) {
  
  diag_data <- patient_f %>% filter(Diagnosis == dx)
  n_sites <- if ("SiteID" %in% names(diag_data)) length(unique(diag_data$SiteID[!is.na(diag_data$SiteID)])) else 0
  
  cat("\n----------------------------------------\n")
  cat(sprintf("Diagnosis: %s | N=%d | Sites=%d\n", dx, nrow(diag_data), n_sites))
  
  t_res <- t.test(diag_data$z_score, mu = 0)
  cat(sprintf("One-sample t-test: t(%.1f)=%.3f, p=%.4g\n", t_res$parameter, t_res$statistic, t_res$p.value))
  
  base_terms <- c("1")
  if ("Age" %in% names(diag_data)) base_terms <- c(base_terms, "Age")
  if ("Sex" %in% names(diag_data)) base_terms <- c(base_terms, "Sex")
  
  fixed_part <- paste(base_terms, collapse = " + ")
  use_mixed  <- ("SiteID" %in% names(diag_data)) && (n_sites > 1)
  
  model_type <- NA
  model_intercept <- NA
  model_se <- NA
  model_t <- NA
  model_p <- NA
  age_beta <- NA; age_p <- NA
  sex_beta <- NA; sex_p <- NA
  
  if (use_mixed) {
    model_type <- "Mixed model: z ~ Age + Sex + (1|SiteID)"
    fml <- as.formula(paste0("z_score ~ ", fixed_part, " + (1|SiteID)"))
    
    tryCatch({
      m <- lmer(fml, data = diag_data)
      ct <- coef(summary(m))
      
      model_intercept <- ct["(Intercept)", "Estimate"]
      model_se        <- ct["(Intercept)", "Std. Error"]
      model_t         <- ct["(Intercept)", "t value"]
      model_p         <- ct["(Intercept)", "Pr(>|t|)"]
      
      if ("Age" %in% colnames(diag_data) && "Age" %in% rownames(ct)) {
        age_beta <- ct["Age", "Estimate"]; age_p <- ct["Age", "Pr(>|t|)"]
      }
      sex_rows <- rownames(ct)[grepl("^Sex", rownames(ct))]
      if (length(sex_rows) > 0) {
        sex_beta <- ct[sex_rows[1], "Estimate"]; sex_p <- ct[sex_rows[1], "Pr(>|t|)"]
      }
      
      cat(sprintf("Model intercept (mean shift): %.3f (SE=%.3f), p=%.4g\n",
                  model_intercept, model_se, model_p))
    }, error = function(e){
      cat(sprintf("Mixed model failed: %s\nFalling back to lm.\n", e$message))
      use_mixed <<- FALSE
    })
  }
  
  if (!use_mixed) {
    model_type <- "Linear regression: z ~ Age + Sex"
    fml <- as.formula(paste0("z_score ~ ", fixed_part))
    
    m <- lm(fml, data = diag_data)
    ct <- coef(summary(m))
    
    model_intercept <- ct["(Intercept)", "Estimate"]
    model_se        <- ct["(Intercept)", "Std. Error"]
    model_t         <- ct["(Intercept)", "t value"]
    model_p         <- ct["(Intercept)", "Pr(>|t|)"]
    
    if ("Age" %in% rownames(ct)) {
      age_beta <- ct["Age", "Estimate"]; age_p <- ct["Age", "Pr(>|t|)"]
    }
    sex_rows <- rownames(ct)[grepl("^Sex", rownames(ct))]
    if (length(sex_rows) > 0) {
      sex_beta <- ct[sex_rows[1], "Estimate"]; sex_p <- ct[sex_rows[1], "Pr(>|t|)"]
    }
    
    cat(sprintf("Model intercept (mean shift): %.3f (SE=%.3f), p=%.4g\n",
                model_intercept, model_se, model_p))
  }
  
  results_list[[as.character(dx)]] <- data.frame(
    Diagnosis = as.character(dx),
    N = nrow(diag_data),
    N_sites = ifelse("SiteID" %in% names(diag_data), n_sites, NA),
    Mean_Z_simple = mean(diag_data$z_score, na.rm=TRUE),
    SD_Z = sd(diag_data$z_score, na.rm=TRUE),
    T_test_t = unname(t_res$statistic),
    T_test_p = t_res$p.value,
    Model_type = model_type,
    Model_intercept = model_intercept,
    Model_SE = model_se,
    Model_t  = model_t,
    Model_p  = model_p,
    Age_effect = age_beta,
    Age_p = age_p,
    Sex_effect = sex_beta,
    Sex_p = sex_p,
    Abnormal_percent = mean(abs(diag_data$z_score) > Zthr, na.rm=TRUE) * 100
  )
}

results_df <- do.call(rbind, results_list) %>% as.data.frame()
rownames(results_df) <- NULL

cat("\n=== Summary of model results ===\n")
print(results_df[, c("Diagnosis","N","Model_type","Model_intercept","Model_SE","Model_p","Abnormal_percent")])

out_res <- file.path(OutDir, "disease_effect_results.csv")
write.csv(results_df, out_res, row.names = FALSE)
cat(sprintf("\nSaved results table to: %s\n", out_res))

dx_order <- desc_stats$Diagnosis
patient_f$Diagnosis <- factor(patient_f$Diagnosis, levels = dx_order)

patient_f <- patient_f %>%
  mutate(
    Abnormal = ifelse(abs(z_score) > Zthr, "Abnormal", "Normal")
  )

sig_df <- results_df %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = levels(patient_f$Diagnosis)),
    p_use = T_test_p,
    sig_lab = case_when(
      is.na(p_use)          ~ "",
      p_use < 0.001         ~ "***",
      p_use < 0.01          ~ "**",
      p_use < 0.05          ~ "*",
      TRUE                  ~ ""
    )
  ) %>%
  filter(sig_lab != "")

y_pos_df <- patient_f %>%
  group_by(Diagnosis) %>%
  summarise(
    y_max = max(z_score, na.rm = TRUE),
    y_min = min(z_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_range = y_max - y_min,
    y_range = ifelse(y_range <= 0 | is.na(y_range), 1, y_range),
    y_pos = y_max + 0.10 * y_range
  )

sig_df <- sig_df %>%
  left_join(y_pos_df, by = "Diagnosis")

dx_levels <- levels(patient_f$Diagnosis)

base_colors <- c(
  "mediumpurple3","plum3","orchid3","skyblue2","steelblue3",
  "lightcoral","salmon3","tan3","darkseagreen3","goldenrod3"
)

if (length(dx_levels) <= length(base_colors)) {
  disease_colors <- setNames(base_colors[seq_along(dx_levels)], dx_levels)
} else {
  disease_colors <- setNames(colorRampPalette(base_colors)(length(dx_levels)), dx_levels)
}

RefLineColor <- "gray55"

y_bottom <- min(patient_f$z_score, na.rm = TRUE)

if (nrow(sig_df) > 0) {
  y_top <- max(sig_df$y_pos, na.rm = TRUE)
} else {
  y_top <- max(patient_f$z_score, na.rm = TRUE)
}

y_pad <- 0.05 * (y_top - y_bottom)
if (!is.finite(y_pad) || y_pad <= 0) y_pad <- 0.5

y_lim <- c(y_bottom - y_pad, y_top + y_pad)

p_rain <- ggplot(
  patient_f,
  aes(x = Diagnosis, y = z_score, fill = Diagnosis)
) +
  ggrain::geom_rain(
    jitter_side = "right",
    point_size = PointSize,
    alpha = 0.65,
    violin.args = list(alpha = 0.45, linewidth = LineWidth),
    boxplot.args = list(width = 0.15, outlier.shape = NA)
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "solid",
    linewidth = LineWidth,
    color = RefLineColor
  ) +
  geom_hline(
    yintercept = c(-Zthr, Zthr),
    linetype = "dashed",
    linewidth = LineWidth,
    color = RefLineColor
  ) +
  geom_text(
    data = sig_df,
    aes(x = Diagnosis, y = y_pos, label = sig_lab),
    inherit.aes = FALSE,
    family = BaseFamily,
    size = 8
  ) +
  scale_fill_manual(values = disease_colors) +
  coord_cartesian(ylim = y_lim, clip = "off") +
  labs(
    title = paste0(PhenotypeName, " z-scores by diagnosis"),
    x = "Diagnosis",
    y = "z-score"
  ) +
  theme_classic(base_family = BaseFamily) +
  theme(
    plot.title = element_text(size = TitleCex, hjust = 0.5),
    axis.title = element_text(size = AxisCex),
    axis.text = element_text(size = AxisCex - 1),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    plot.margin = margin(10, 20, 10, 10)
  )

print(p_rain)

ggsave(
  filename = file.path(OutDir, "raincloud_z_by_diagnosis.png"),
  plot = p_rain,
  width = FigW,
  height = FigH,
  dpi = DPI,
  bg = "white"
)

cat("\nRaincloud plot with significance stars saved to OutDir.\n")
