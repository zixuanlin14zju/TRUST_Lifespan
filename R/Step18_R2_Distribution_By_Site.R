# Step18_R2_Distribution_By_Site.R
# Public analysis entry point. Start from the repository root, or use:
#   Rscript run_step.R 18
source("R/lib/oef_utils.R")
oef_start_step("Step18")

library(readxl)

DataPath <- oef_input("trust_hc_hct_ya_sensitivity.xlsx")
OutDir  <- oef_output("r2_by_site")
dir.create(OutDir, recursive = TRUE, showWarnings = FALSE)

M_HC <- read_excel(DataPath, sheet = "Sheet1")

# Resting venous blood R2, assuming T2 is in ms
oef_require_columns(M_HC, c("SiteID", "T2"))
if (any(!is.na(M_HC$T2) & (!is.finite(M_HC$T2) | M_HC$T2 <= 0))) stop("T2 must be positive milliseconds.")
M_HC$R2 <- 1000 / M_HC$T2

# Plot settings: identical to original Figure 1
FontFamily <- "sans"
TextCex    <- 1.8
AxisCex    <- 1.4
TitleCex   <- 1.5
LabCex     <- 1.5
DPI        <- 600
FigW_in    <- 10
FigH_in    <- 7

open_png <- function(filename, w = FigW_in, h = FigH_in, res = DPI){
  f <- file.path(OutDir, filename)
  ok <- FALSE
  
  try({
    png(f, width = w, height = h, units = "in",
        res = res, type = "cairo")
    ok <- TRUE
  }, silent = TRUE)
  
  if (!ok)
    png(f, width = w, height = h, units = "in", res = res)
  
  par(
    family = FontFamily,
    cex    = TextCex,
    mar    = c(5, 5, 4, 2) + 0.2,
    mgp    = c(2.6, 0.9, 0),
    tcl    = -0.3,
    bty    = "l"
  )
}

close_png <- function(){
  try(dev.off(), silent = TRUE)
}

plot_labs <- function(xlab, ylab, main){
  title(main = main, cex.main = TitleCex)
  title(xlab = xlab, ylab = ylab, cex.lab = LabCex)
  par(cex.axis = AxisCex)
}

# Same site order as the original Figure 1
site_levels <- sort(unique(M_HC$SiteID))
M_HC$SiteID_plot <- factor(M_HC$SiteID, levels = site_levels)

n_sites <- length(site_levels)

# Warm-color gradient for R2
r2_cols <- colorRampPalette(
  c("#fee6ce", "#a63603")
)(n_sites)

open_png("01_05_R2_by_site.png")

boxplot(
  R2 ~ SiteID_plot,
  data = M_HC,
  col = r2_cols,
  border = "black",
  xaxt = "n",
  cex = 0.6,
  xlab = "",
  ylab = "",
  main = ""
)

axis(
  1,
  at = seq_along(site_levels),
  labels = site_levels,
  las = 2,
  cex.axis = AxisCex * 0.7
)

plot_labs(
  "Site ID",
  expression(R[2]~(s^{-1})),
  "Resting venous blood R2 by site"
)

close_png()
