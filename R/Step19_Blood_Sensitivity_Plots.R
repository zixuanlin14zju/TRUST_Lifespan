# Step19_Blood_Sensitivity_Plots.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 19
source("R/lib/oef_utils.R")
oef_start_step("Step19")

# ============================================================
# Sex-specific plotting for Hct / Ya sensitivity analyses
#
# This script reads the existing GAMLSS output tables and generates plots for:
#   1) median OEF trajectory
#   2) growth-rate (first derivative) trajectory
#
# For BOTH sensitivity analyses:
#   - Hct sensitivity
#   - Ya sensitivity
#
# Each figure contains three panels:
#   Overall / Female / Male
# ============================================================

# ============================================================
# 00) Packages
# ============================================================

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "scales"
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
  library(readr)
  library(stringr)
  library(scales)
})

# ============================================================
# 01) User settings
# ============================================================

# --- Change these paths if needed ---
HctInput <- oef_output("GAMLSS_HctSensitivity")
YaInput  <- oef_output("GAMLSS_YaSensitivity_Hybrid_FINAL")

OutDir <- oef_output("blood_sensitivity_plots")
dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)

# Export settings
DPI <- 600
FigWidth <- 10.5
FigHeight <- 11.0

# Keep the same y-axis across panels within the same figure
# (e.g., Overall / Female / Male share one common y-range).
SameYAxisWithinFigure <- TRUE

# If you only want ages >= 1 year shown in the main panels, set to 1.
# Keep at 0.25 to show the neonatal range.
PlotAgeMin <- 0.25

# Optional zoomed adult views
MakeAdultOnlyPlots <- TRUE
AdultAgeMin <- 18

# ============================================================
# 02) Helpers
# ============================================================

is_zip_path <- function(x) {
  grepl("\\.zip$", x, ignore.case = TRUE)
}

list_zip_files <- function(zip_path) {
  utils::unzip(zip_path, list = TRUE)$Name
}

read_csv_from_folder_or_zip <- function(input_path, candidate_files) {

  # candidate_files: ordered vector of possible inner/relative file names
  # returns the first match found

  if (!is_zip_path(input_path)) {

    for (f in candidate_files) {
      full <- file.path(input_path, f)
      if (file.exists(full)) {
        return(readr::read_csv(full, show_col_types = FALSE))
      }
    }

    stop(
      "Could not find any of these files in folder:\n",
      paste(candidate_files, collapse = "\n"),
      "\nFolder: ", input_path
    )

  } else {

    zfiles <- list_zip_files(input_path)

    for (f in candidate_files) {
      hit <- zfiles[zfiles == f]
      if (length(hit) == 1) {
        return(readr::read_csv(unz(input_path, hit), show_col_types = FALSE))
      }
    }

    stop(
      "Could not find any of these files in zip:\n",
      paste(candidate_files, collapse = "\n"),
      "\nZip: ", input_path
    )
  }
}

save_plot <- function(p, filename, width = FigWidth, height = FigHeight) {
  ggplot2::ggsave(
    filename = file.path(OutDir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = DPI,
    bg = "white"
  )
}

theme_trust <- function() {
  theme_bw(base_size = 15) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold", size = 20),
      legend.title = element_blank(),
      legend.position = "top",
      legend.box = "horizontal",
      legend.key.width = unit(2.0, "lines"),
      legend.text = element_text(size = 15),
      axis.title = element_text(face = "bold", size = 22),
      plot.title = element_text(face = "bold", size = 20, hjust = 0),
      plot.subtitle = element_text(size = 18, hjust = 0),
      axis.text = element_text(color = "black", size = 20)
    )
}

scenario_factor <- function(x, levels_vec) {
  factor(x, levels = levels_vec)
}

pivot_curve_data <- function(df) {

  # expects columns:
  # Scenario, Age,
  # Median_Overall, Median_Female, Median_Male,
  # Growth_Overall, Growth_Female, Growth_Male (for Ya)
  # OR trajectory/growth already separate

  out <- df %>%
    select(any_of(c(
      "Scenario", "Age",
      "Median_Overall", "Median_Female", "Median_Male",
      "Growth_Overall", "Growth_Female", "Growth_Male"
    )))

  traj_long <- out %>%
    select(Scenario, Age, starts_with("Median_")) %>%
    pivot_longer(
      cols = starts_with("Median_"),
      names_to = "Panel",
      values_to = "Value"
    ) %>%
    mutate(
      CurveType = "Trajectory",
      Panel = recode(
        Panel,
        Median_Overall = "Overall",
        Median_Female  = "Female",
        Median_Male    = "Male"
      )
    )

  growth_long <- out %>%
    select(Scenario, Age, starts_with("Growth_")) %>%
    pivot_longer(
      cols = starts_with("Growth_"),
      names_to = "Panel",
      values_to = "Value"
    ) %>%
    mutate(
      CurveType = "GrowthRate",
      Panel = recode(
        Panel,
        Growth_Overall = "Overall",
        Growth_Female  = "Female",
        Growth_Male    = "Male"
      )
    )

  list(
    trajectory = traj_long,
    growth = growth_long
  )
}

build_threepanel_plot <- function(
  dat,
  palette,
  linetypes,
  title,
  subtitle,
  ylab,
  filename,
  x_min = PlotAgeMin,
  x_max = NULL,
  add_zero_line = FALSE,
  same_y_across_panels = TRUE
) {

  dat <- dat %>%
    filter(Age >= x_min)

  if (!is.null(x_max)) {
    dat <- dat %>% filter(Age <= x_max)
  }

  facet_scales <- if (same_y_across_panels) "fixed" else "free_y"

  p <- ggplot(
    dat,
    aes(x = Age, y = Value, color = Scenario, linetype = Scenario)
  ) +
    geom_line(linewidth = 1.3) +
    facet_wrap(~Panel, ncol = 1, scales = facet_scales) +
    scale_color_manual(
      values = palette,
      labels = c(
        "Original" = "Original",
        "RegionMatched" = "Region-matched",
        "Zierk" = "Zierk",
        "MahlknechtMean" = "Mahlknecht",
        "CohortEmpirical" = "Cohort-empirical",
        "LuHybrid" = "Lu",
        "CrapoHybrid" = "Crapo"
      ),
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = linetypes,
      labels = c(
        "Original" = "Original",
        "RegionMatched" = "Region-matched",
        "Zierk" = "Zierk",
        "MahlknechtMean" = "Mahlknecht",
        "CohortEmpirical" = "Cohort-empirical",
        "LuHybrid" = "Lu",
        "CrapoHybrid" = "Crapo"
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 8),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Age (years)",
      y = ylab
    ) +
    theme_trust()

  if (add_zero_line) {
    p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey50")
  }

  save_plot(p, filename)
  p
}

# ============================================================
# 03) Read Hct sensitivity outputs
# ============================================================

HctTraj <- read_csv_from_folder_or_zip(
  HctInput,
  c(
    "trajectory_by_age.csv",
    "GAMLSS_HctSensitivity/trajectory_by_age.csv"
  )
)

HctGrowth <- read_csv_from_folder_or_zip(
  HctInput,
  c(
    "growth_rate_by_age.csv",
    "GAMLSS_HctSensitivity/growth_rate_by_age.csv"
  )
)

Hct <- HctTraj %>%
  left_join(
    HctGrowth,
    by = c("Scenario", "Age")
  )

hct_levels <- c(
  "Original",
  "RegionMatched",
  "Zierk",
  "MahlknechtMean",
  "CohortEmpirical"
)

Hct$Scenario <- scenario_factor(Hct$Scenario, hct_levels)

HctLong <- pivot_curve_data(Hct)

hct_palette <- c(
  Original = "black",
  RegionMatched = "#1b9e77",
  Zierk = "#377eb8",
  MahlknechtMean = "#d95f02",
  CohortEmpirical = "#7570b3"
)

hct_linetypes <- c(
  Original = "solid",
  RegionMatched = "solid",
  Zierk = "dashed",
  MahlknechtMean = "dotdash",
  CohortEmpirical = "longdash"
)

# ============================================================
# 04) Read Ya sensitivity outputs
# ============================================================

Ya <- read_csv_from_folder_or_zip(
  YaInput,
  c(
    "trajectory_and_growth_by_age.csv",
    "01_FullLifespan_Original_vs_LuHybrid_vs_CrapoHybrid/trajectory_and_growth_by_age.csv"
  )
)

ya_levels <- c(
  "Original",
  "LuHybrid",
  "CrapoHybrid"
)

Ya$Scenario <- scenario_factor(Ya$Scenario, ya_levels)

YaLong <- pivot_curve_data(Ya)

ya_palette <- c(
  Original = "black",
  LuHybrid = "#1b9e77",
  CrapoHybrid = "#d95f02"
)

ya_linetypes <- c(
  Original = "solid",
  LuHybrid = "solid",
  CrapoHybrid = "dashed"
)

# ============================================================
# 05) Main full-lifespan plots
# ============================================================

# ----- Hct trajectory -----
build_threepanel_plot(
  dat = HctLong$trajectory,
  palette = hct_palette,
  linetypes = hct_linetypes,
  title = "Hematocrit sensitivity: fitted OEF trajectory",
  subtitle = "Panels show overall, female, and male normative trajectories",
  ylab = "Median OEF",
  filename = "HCT_01_Trajectory_Overall_Female_Male.png",
  same_y_across_panels = SameYAxisWithinFigure
)

# ----- Hct growth-rate -----
build_threepanel_plot(
  dat = HctLong$growth,
  palette = hct_palette,
  linetypes = hct_linetypes,
  title = "Hematocrit sensitivity: OEF growth-rate curve",
  subtitle = "First derivative of the fitted normative OEF trajectory",
  ylab = "dOEF / dAge",
  filename = "HCT_02_GrowthRate_Overall_Female_Male.png",
  add_zero_line = TRUE,
  same_y_across_panels = SameYAxisWithinFigure
)

# ----- Ya trajectory -----
build_threepanel_plot(
  dat = YaLong$trajectory,
  palette = ya_palette,
  linetypes = ya_linetypes,
  title = "Arterial oxygen saturation sensitivity: fitted OEF trajectory",
  subtitle = "Panels show overall, female, and male normative trajectories",
  ylab = "Median OEF",
  filename = "YA_01_Trajectory_Overall_Female_Male.png",
  same_y_across_panels = SameYAxisWithinFigure
)

# ----- Ya growth-rate -----
build_threepanel_plot(
  dat = YaLong$growth,
  palette = ya_palette,
  linetypes = ya_linetypes,
  title = "Arterial oxygen saturation sensitivity: OEF growth-rate curve",
  subtitle = "First derivative of the fitted normative OEF trajectory",
  ylab = "dOEF / dAge",
  filename = "YA_02_GrowthRate_Overall_Female_Male.png",
  add_zero_line = TRUE,
  same_y_across_panels = SameYAxisWithinFigure
)

# ============================================================
# 06) Optional adult-only versions
# ============================================================

if (MakeAdultOnlyPlots) {

  build_threepanel_plot(
    dat = HctLong$trajectory,
    palette = hct_palette,
    linetypes = hct_linetypes,
    title = "Hematocrit sensitivity: adult OEF trajectory",
    subtitle = paste0("Restricted to age ≥ ", AdultAgeMin, " years"),
    ylab = "Median OEF",
    filename = "HCT_03_AdultOnly_Trajectory_Overall_Female_Male.png",
    x_min = AdultAgeMin,
    same_y_across_panels = SameYAxisWithinFigure
  )

  build_threepanel_plot(
    dat = HctLong$growth,
    palette = hct_palette,
    linetypes = hct_linetypes,
    title = "Hematocrit sensitivity: adult OEF growth-rate curve",
    subtitle = paste0("Restricted to age ≥ ", AdultAgeMin, " years"),
    ylab = "dOEF / dAge",
    filename = "HCT_04_AdultOnly_GrowthRate_Overall_Female_Male.png",
    x_min = AdultAgeMin,
    add_zero_line = TRUE,
    same_y_across_panels = SameYAxisWithinFigure
  )

  build_threepanel_plot(
    dat = YaLong$trajectory,
    palette = ya_palette,
    linetypes = ya_linetypes,
    title = "Arterial oxygen saturation sensitivity: adult OEF trajectory",
    subtitle = paste0("Restricted to age ≥ ", AdultAgeMin, " years"),
    ylab = "Median OEF",
    filename = "YA_03_AdultOnly_Trajectory_Overall_Female_Male.png",
    x_min = AdultAgeMin,
    same_y_across_panels = SameYAxisWithinFigure
  )

  build_threepanel_plot(
    dat = YaLong$growth,
    palette = ya_palette,
    linetypes = ya_linetypes,
    title = "Arterial oxygen saturation sensitivity: adult OEF growth-rate curve",
    subtitle = paste0("Restricted to age ≥ ", AdultAgeMin, " years"),
    ylab = "dOEF / dAge",
    filename = "YA_04_AdultOnly_GrowthRate_Overall_Female_Male.png",
    x_min = AdultAgeMin,
    add_zero_line = TRUE,
    same_y_across_panels = SameYAxisWithinFigure
  )
}

# ============================================================
# 07) Also make a compact 2-panel sex-only version
#     (sometimes cleaner for Supplement figures)
# ============================================================

build_two_sex_plot <- function(
  dat,
  palette,
  linetypes,
  title,
  subtitle,
  ylab,
  filename,
  x_min = PlotAgeMin,
  add_zero_line = FALSE,
  same_y_across_panels = TRUE
) {

  dat2 <- dat %>%
    filter(Panel %in% c("Female", "Male")) %>%
    filter(Age >= x_min)

  facet_scales <- if (same_y_across_panels) "fixed" else "free_y"

  p <- ggplot(
    dat2,
    aes(x = Age, y = Value, color = Scenario, linetype = Scenario)
  ) +
    geom_line(linewidth = 1.3) +
    facet_wrap(~Panel, ncol = 1, scales = facet_scales) +
    scale_color_manual(
      values = palette,
      labels = c(
        "Original" = "Original",
        "RegionMatched" = "Region-matched",
        "Zierk" = "Zierk",
        "MahlknechtMean" = "Mahlknecht",
        "CohortEmpirical" = "Cohort-empirical",
        "LuHybrid" = "Lu",
        "CrapoHybrid" = "Crapo"
      ),
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = linetypes,
      labels = c(
        "Original" = "Original",
        "RegionMatched" = "Region-matched",
        "Zierk" = "Zierk",
        "MahlknechtMean" = "Mahlknecht",
        "CohortEmpirical" = "Cohort-empirical",
        "LuHybrid" = "Lu",
        "CrapoHybrid" = "Crapo"
      ),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = pretty_breaks(n = 8),
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Age (years)",
      y = ylab
    ) +
    theme_trust()

  if (add_zero_line) {
    p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey50")
  }

  save_plot(p, filename, width = 9.0, height = 8.2)
  p
}

build_two_sex_plot(
  dat = HctLong$trajectory,
  palette = hct_palette,
  linetypes = hct_linetypes,
  title = "Hematocrit sensitivity: female vs male OEF trajectory",
  subtitle = "Sex-specific normative curves only",
  ylab = "Median OEF",
  filename = "HCT_05_SexOnly_Trajectory.png",
  same_y_across_panels = SameYAxisWithinFigure
)

build_two_sex_plot(
  dat = HctLong$growth,
  palette = hct_palette,
  linetypes = hct_linetypes,
  title = "Hematocrit sensitivity: female vs male OEF growth rate",
  subtitle = "Sex-specific first-derivative curves only",
  ylab = "dOEF / dAge",
  filename = "HCT_06_SexOnly_GrowthRate.png",
  add_zero_line = TRUE,
  same_y_across_panels = SameYAxisWithinFigure
)

build_two_sex_plot(
  dat = YaLong$trajectory,
  palette = ya_palette,
  linetypes = ya_linetypes,
  title = "Arterial oxygen saturation sensitivity: female vs male OEF trajectory",
  subtitle = "Sex-specific normative curves only",
  ylab = "Median OEF",
  filename = "YA_05_SexOnly_Trajectory.png",
  same_y_across_panels = SameYAxisWithinFigure
)

build_two_sex_plot(
  dat = YaLong$growth,
  palette = ya_palette,
  linetypes = ya_linetypes,
  title = "Arterial oxygen saturation sensitivity: female vs male OEF growth rate",
  subtitle = "Sex-specific first-derivative curves only",
  ylab = "dOEF / dAge",
  filename = "YA_06_SexOnly_GrowthRate.png",
  add_zero_line = TRUE,
  same_y_across_panels = SameYAxisWithinFigure
)

cat("\n============================================================\n")
cat("Pretty plotting completed.\n")
cat("Output directory:\n")
cat(OutDir, "\n")
cat("============================================================\n")
