# Step20_Blood_Assumption_Plots.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 20
source("R/lib/oef_utils.R")
oef_start_step("Step20")

# ============================================================
# Supplementary Figure: Hct- and Ya-age assumptions
#
# Purpose:
#   Visualize the blood-parameter assumptions used in the
#   Hct and Ya sensitivity analyses.
#
# Figure layout:
#   A. Hct vs age - Male
#   B. Hct vs age - Female
#   C. Ya vs age
#
# Important interpretation:
#   1) Hct curves show the MODEL-BASED/reference assumptions.
#      Individual measured-Hct overrides are not plotted.
#   2) The Region-matched Hct curve shown here is for Chinese sites.
#      BIOCARD (Site 11) uses the Zierk Hct reference.
#   3) For Ya, participants <18 y retained measured Ya whenever
#      available. The Lu/Crapo equations were used only as pediatric
#      fallbacks, but were applied to all adults (>=18 y).
# ============================================================


# ============================================================
# 00) Packages
# ============================================================

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "scales",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Please install package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nExample:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})


# ============================================================
# 01) User settings
# ============================================================

OutDir <- oef_output("blood_sensitivity_plots")
dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)

DPI <- 600

# Start at 0.25 y because the PRINCE continuous reference used
# in the Region-matched model begins at 3 months.
AgeMin <- 0.25
AgeMax <- 93

# Fine age grid for smooth interpolation curves
AgeStep <- 0.01
age_grid <- seq(AgeMin, AgeMax, by = AgeStep)

FigWidth  <- 10.5
FigHeight <- 12.0

# Larger font because the three panels will be combined
BaseFontSize <- 15


# ============================================================
# 02) Plot theme
# ============================================================

theme_assumption <- function() {
  theme_bw(base_size = BaseFontSize) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title = element_text(face = "bold", size = BaseFontSize + 1),
      axis.text = element_text(size = BaseFontSize - 1, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = BaseFontSize - 1),
      legend.key.width = unit(2.0, "lines"),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = BaseFontSize + 2),
      plot.subtitle = element_text(size = BaseFontSize - 1),
      plot.tag = element_text(face = "bold", size = BaseFontSize + 4),
      plot.margin = margin(8, 12, 8, 8)
    )
}


# ============================================================
# 03) Hct model definitions
# ============================================================

# ------------------------------------------------------------
# 03a) PRINCE pediatric midpoint proxy used for Chinese sites
#      in the Region-matched model
#
# Each anchor is:
#   (continuous lower RI + continuous upper RI) / 2
# from PRINCE Supplementary Table S19.
# Participant-specific values were linearly interpolated by age.
# ------------------------------------------------------------

prince_age <- c(
  0.25, 0.5,
  1, 2, 3, 4, 5, 6, 7, 8, 9,
  10, 11, 12, 13, 14, 15, 16, 17, 18, 19
)

prince_male <- c(
  35.0, 35.0,
  36.0, 37.0, 38.0, 38.5, 38.5, 39.0, 39.5, 39.5, 40.0,
  40.5, 41.0, 42.0, 43.0, 44.5, 45.0, 46.0, 46.0, 46.0, 46.0
)

prince_female <- c(
  35.5, 35.5,
  36.5, 37.5, 37.5, 38.5, 38.5, 39.0, 39.5, 39.5, 40.5,
  40.5, 41.0, 41.0, 41.0, 41.5, 41.5, 41.0, 41.0, 41.0, 40.5
)


# ------------------------------------------------------------
# 03b) Zierk/PEDREF pediatric digitized 50th-percentile anchors
# ------------------------------------------------------------

zierk_age <- c(
  0,
  30 / 365.25,
  60 / 365.25,
  90 / 365.25,
  120 / 365.25,
  180 / 365.25,
  1, 2, 4, 6, 8, 10, 12, 14, 16, 18
)

zierk_male <- c(
  53.5, 40.5, 31.5, 31.5, 34.0, 34.0,
  34.5, 35.0, 36.0, 37.0, 38.0, 38.5, 40.0, 41.0, 42.5, 45.0
)

zierk_female <- c(
  53.5, 40.5, 31.5, 31.5, 34.0, 34.0,
  34.5, 35.0, 36.0, 37.0, 38.0, 38.0, 38.0, 38.0, 38.0, 38.0
)


# ------------------------------------------------------------
# 03c) Helper functions
# ------------------------------------------------------------

linear_bridge <- function(age, age0, age1, value0, value1) {
  value0 + (age - age0) / (age1 - age0) * (value1 - value0)
}


region_matched_hct <- function(age, sex = c("Male", "Female")) {

  sex <- match.arg(sex)

  out <- rep(NA_real_, length(age))

  # PRINCE: 0.25-19 y, exact-age linear interpolation
  idx <- age >= 0.25 & age <= 19
  if (any(idx)) {
    vals <- if (sex == "Male") prince_male else prince_female
    out[idx] <- approx(
      x = prince_age,
      y = vals,
      xout = age[idx],
      method = "linear",
      rule = 1
    )$y
  }

  # 19-20 y bridge from PRINCE age-19 value to Wu 20-29 y mean
  idx <- age > 19 & age < 20
  if (any(idx)) {
    if (sex == "Male") {
      out[idx] <- linear_bridge(age[idx], 19, 20, 46.0, 46.2)
    } else {
      out[idx] <- linear_bridge(age[idx], 19, 20, 40.5, 39.7)
    }
  }

  # Wu adult age-group means
  if (sex == "Male") {
    out[age >= 20 & age < 30] <- 46.2
    out[age >= 30 & age < 40] <- 46.2
    out[age >= 40 & age < 50] <- 46.2
    out[age >= 50 & age < 60] <- 45.5
    out[age >= 60 & age < 70] <- 44.8
    out[age >= 70 & age < 80] <- 44.1

    # Chinese older-adult reference
    out[age >= 80] <- 43.5

  } else {

    out[age >= 20 & age < 30] <- 39.7
    out[age >= 30 & age < 40] <- 40.1
    out[age >= 40 & age < 50] <- 40.1
    out[age >= 50 & age < 60] <- 40.7
    out[age >= 60 & age < 70] <- 40.7
    out[age >= 70 & age < 80] <- 40.7

    # Chinese older-adult reference
    out[age >= 80] <- 39.8
  }

  out
}


zierk_hct <- function(age, sex = c("Male", "Female")) {

  sex <- match.arg(sex)
  out <- rep(NA_real_, length(age))

  # Pediatric curve: digitized 50th percentile, interpolated by exact age
  idx <- age <= 18
  if (any(idx)) {
    vals <- if (sex == "Male") zierk_male else zierk_female
    out[idx] <- approx(
      x = zierk_age,
      y = vals,
      xout = age[idx],
      method = "linear",
      rule = 1
    )$y
  }

  # 18-20 y bridge to adult Zierk medians
  idx <- age > 18 & age < 20
  if (any(idx)) {
    if (sex == "Male") {
      out[idx] <- linear_bridge(age[idx], 18, 20, 45.0, 44.5)
    } else {
      out[idx] <- linear_bridge(age[idx], 18, 20, 38.0, 39.5)
    }
  }

  # Adult age-specific medians
  if (sex == "Male") {
    out[age >= 20 & age < 30] <- 44.5
    out[age >= 30 & age < 40] <- 44.3
    out[age >= 40 & age < 50] <- 44.1
    out[age >= 50 & age < 60] <- 43.8
    out[age >= 60 & age < 70] <- 43.4
    out[age >= 70 & age < 80] <- 42.3
    out[age >= 80 & age < 90] <- 41.4
    out[age >= 90] <- 39.7
  } else {
    out[age >= 20 & age < 30] <- 39.5
    out[age >= 30 & age < 40] <- 39.5
    out[age >= 40 & age < 50] <- 39.8
    out[age >= 50 & age < 60] <- 40.3
    out[age >= 60 & age < 70] <- 40.5
    out[age >= 70 & age < 80] <- 40.0
    out[age >= 80 & age < 90] <- 39.4
    out[age >= 90] <- 38.8
  }

  out
}


mahlknecht_hct <- function(age, sex = c("Male", "Female")) {

  sex <- match.arg(sex)
  out <- rep(NA_real_, length(age))

  # Pediatric fill in the FINAL analysis: same as Zierk
  idx <- age < 18
  out[idx] <- zierk_hct(age[idx], sex)

  # Mahlknecht reconstructed means from age 18 onward
  if (sex == "Male") {
    out[age >= 18 & age < 20] <- 44.5
    out[age >= 20 & age < 30] <- 44.0
    out[age >= 30 & age < 40] <- 43.5
    out[age >= 40 & age < 50] <- 42.5
    out[age >= 50 & age < 60] <- 41.0
    out[age >= 60 & age < 70] <- 41.0
    out[age >= 70 & age < 80] <- 37.5
    out[age >= 80 & age < 90] <- 35.5
    out[age >= 90] <- 36.0
  } else {
    out[age >= 18 & age < 20] <- 39.0
    out[age >= 20 & age < 30] <- 39.5
    out[age >= 30 & age < 40] <- 38.5
    out[age >= 40 & age < 50] <- 38.5
    out[age >= 50 & age < 60] <- 40.0
    out[age >= 60 & age < 70] <- 38.0
    out[age >= 70 & age < 80] <- 36.5
    out[age >= 80 & age < 90] <- 35.5
    out[age >= 90] <- 33.5
  }

  out
}


cohort_empirical_hct <- function(age, sex = c("Male", "Female")) {

  sex <- match.arg(sex)
  out <- rep(NA_real_, length(age))

  pediatric_max <- 12.616
  adult_min <- 30

  # Pediatric and adult regressions were fitted separately;
  # within each regression, age and sex were predictors.
  if (sex == "Male") {

    pediatric_fun <- function(a) {
      39.52999 + 0.08404 * a
    }

    adult_fun <- function(a) {
      47.12804 - 0.08514 * a
    }

  } else {

    pediatric_fun <- function(a) {
      39.69531 + 0.08404 * a
    }

    adult_fun <- function(a) {
      44.07132 - 0.08514 * a
    }
  }

  # Within observed pediatric measured-data range
  idx <- age <= pediatric_max
  out[idx] <- pediatric_fun(age[idx])

  # Bridge across the gap in measured-data coverage
  idx <- age > pediatric_max & age < adult_min
  if (any(idx)) {
    out[idx] <- linear_bridge(
      age = age[idx],
      age0 = pediatric_max,
      age1 = adult_min,
      value0 = pediatric_fun(pediatric_max),
      value1 = adult_fun(adult_min)
    )
  }

  # Within adult measured-data range
  idx <- age >= adult_min
  out[idx] <- adult_fun(age[idx])

  out
}


fixed_hct <- function(age, sex = c("Male", "Female")) {
  sex <- match.arg(sex)
  rep(ifelse(sex == "Male", 42.0, 40.0), length(age))
}


# ============================================================
# 04) Build Hct plotting data
# ============================================================

make_hct_data <- function(sex_label) {

  tibble(
    Age = age_grid,
    `Fixed 42M/40F` = fixed_hct(age_grid, sex_label),
    `Region-matched` = region_matched_hct(age_grid, sex_label),
    `Zierk` = zierk_hct(age_grid, sex_label),
    `Mahlknecht` = mahlknecht_hct(age_grid, sex_label),
    `Cohort-empirical` = cohort_empirical_hct(age_grid, sex_label)
  ) %>%
    pivot_longer(
      cols = -Age,
      names_to = "Model",
      values_to = "Hct"
    ) %>%
    mutate(
      Sex = sex_label,
      Model = factor(
        Model,
        levels = c(
          "Fixed 42M/40F",
          "Region-matched",
          "Zierk",
          "Mahlknecht",
          "Cohort-empirical"
        )
      )
    )
}

HctPlotData <- bind_rows(
  make_hct_data("Male"),
  make_hct_data("Female")
)


# ============================================================
# 05) Ya model definitions
# ============================================================

# Fixed 98% reference used as the standard assumption when
# participant-specific arterial saturation was unavailable.
ya_fixed <- function(age) {
  rep(98.0, length(age))
}

# Lu et al.
ya_lu <- function(age) {
  99.06 - 0.02 * age
}

# Crapo et al.
ya_crapo <- function(age) {
  97.66 - 0.0296 * age
}


YaPlotData <- tibble(
  Age = age_grid,
  `Fixed 98%` = ya_fixed(age_grid),
  `Lu` = ya_lu(age_grid),
  `Crapo` = ya_crapo(age_grid)
) %>%
  pivot_longer(
    cols = -Age,
    names_to = "Model",
    values_to = "Ya"
  ) %>%
  mutate(
    Model = factor(
      Model,
      levels = c(
        "Fixed 98%",
        "Lu",
        "Crapo"
      )
    )
  )


# ============================================================
# 06) Palettes and linetypes
# ============================================================

hct_palette <- c(
  "Fixed 42M/40F" = "black",
  "Region-matched" = "#1b9e77",
  "Zierk" = "#377eb8",
  "Mahlknecht" = "#d95f02",
  "Cohort-empirical" = "#7570b3"
)

hct_linetypes <- c(
  "Fixed 42M/40F" = "solid",
  "Region-matched" = "solid",
  "Zierk" = "dashed",
  "Mahlknecht" = "dotdash",
  "Cohort-empirical" = "longdash"
)

ya_palette <- c(
  "Fixed 98%" = "black",
  "Lu" = "#1b9e77",
  "Crapo" = "#d95f02"
)

ya_linetypes <- c(
  "Fixed 98%" = "solid",
  "Lu" = "solid",
  "Crapo" = "dashed"
)


# ============================================================
# 07) Hct panels
# ============================================================

build_hct_panel <- function(sex_label) {

  dat <- HctPlotData %>%
    filter(Sex == sex_label)

  ggplot(
    dat,
    aes(
      x = Age,
      y = Hct,
      color = Model,
      linetype = Model
    )
  ) +
    geom_line(linewidth = 1.25, na.rm = TRUE) +
    scale_color_manual(values = hct_palette, drop = FALSE) +
    scale_linetype_manual(values = hct_linetypes, drop = FALSE) +
    scale_x_continuous(
      breaks = seq(0, 90, by = 10),
      limits = c(AgeMin, AgeMax),
      expand = expansion(mult = c(0.005, 0.01))
    ) +
    scale_y_continuous(
      breaks = seq(30, 50, by = 5),
      limits = c(30, 48)
    ) +
    labs(
      title = paste0("Hct assumptions: ", tolower(sex_label)),
      x = "Age (years)",
      y = "Hematocrit (%)"
    ) +
    theme_assumption()
}


p_hct_male <- build_hct_panel("Male")
p_hct_female <- build_hct_panel("Female")


# ============================================================
# 08) Ya panel
# ============================================================

p_ya <- ggplot(
  YaPlotData,
  aes(
    x = Age,
    y = Ya,
    color = Model,
    linetype = Model
  )
) +
  geom_line(linewidth = 1.3) +

  # Age 18 marks the transition at which Lu/Crapo equations
  # were applied to every participant rather than only as fallback.
  geom_vline(
    xintercept = 18,
    linetype = "dotted",
    linewidth = 0.8,
    color = "grey40"
  ) +

  annotate(
    "text",
    x = 19.5,
    y = 99.25,
    label = "Age 18",
    hjust = 0,
    size = 4.6
  ) +

  scale_color_manual(values = ya_palette, drop = FALSE) +
  scale_linetype_manual(values = ya_linetypes, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(0, 90, by = 10),
    limits = c(AgeMin, AgeMax),
    expand = expansion(mult = c(0.005, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(95, 99, by = 1),
    limits = c(94.5, 99.5)
  ) +
  labs(
    title = "Arterial oxygen saturation assumptions",
    subtitle = "For age <18 y, measured Ya was retained when available; Lu/Crapo equations served as fallback assumptions",
    x = "Age (years)",
    y = expression(Y[a]~"(%)")
  ) +
  theme_assumption()


# ============================================================
# 09) Combined supplementary figure
# ============================================================

combined_plot <-
  p_hct_male /
  p_hct_female /
  p_ya +
  plot_annotation(
    tag_levels = "A",
    title = "Blood-parameter assumptions used in sensitivity analyses",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = BaseFontSize + 4,
        hjust = 0
      )
    )
  )


# ============================================================
# 10) Save outputs
# ============================================================

ggsave(
  filename = file.path(OutDir, "SuppFig_Hct_Ya_ModelAssumptions.png"),
  plot = combined_plot,
  width = FigWidth,
  height = FigHeight,
  dpi = DPI,
  bg = "white"
)

ggsave(
  filename = file.path(OutDir, "SuppFig_Hct_Ya_ModelAssumptions.pdf"),
  plot = combined_plot,
  width = FigWidth,
  height = FigHeight,
  device = cairo_pdf,
  bg = "white"
)

# Also save individual panels in case the final layout needs adjustment
ggsave(
  filename = file.path(OutDir, "SuppFig_Hct_Assumptions_Male.png"),
  plot = p_hct_male,
  width = 10,
  height = 5.5,
  dpi = DPI,
  bg = "white"
)

ggsave(
  filename = file.path(OutDir, "SuppFig_Hct_Assumptions_Female.png"),
  plot = p_hct_female,
  width = 10,
  height = 5.5,
  dpi = DPI,
  bg = "white"
)

ggsave(
  filename = file.path(OutDir, "SuppFig_Ya_Assumptions.png"),
  plot = p_ya,
  width = 10,
  height = 5.5,
  dpi = DPI,
  bg = "white"
)


# ============================================================
# 11) Console audit
# ============================================================

cat("\n============================================================\n")
cat("Model-assumption plotting completed.\n")
cat("Output directory:\n")
cat(OutDir, "\n\n")
cat("Important interpretation:\n")
cat("  - Hct curves show model/reference assumptions only.\n")
cat("  - Individual measured-Hct overrides are not plotted.\n")
cat("  - Region-matched curve shown = Chinese sites.\n")
cat("  - BIOCARD uses the Zierk Hct reference.\n")
cat("  - For Ya <18 y, measured Ya was retained when available.\n")
cat("  - Lu/Crapo pediatric curves represent fallback assumptions.\n")
cat("============================================================\n")
