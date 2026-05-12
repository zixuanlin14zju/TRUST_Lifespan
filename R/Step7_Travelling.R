library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(irr)

file_path <- "data/trust_traveling_subjects.xlsx"
sheet_name <- 1
outdir <- "outputs/reproducibility_traveling"

if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

dat <- read_excel(file_path, sheet = sheet_name) %>%
  as.data.frame() %>%
  mutate(name = as.character(name)) 

calc_icc <- function(x, y){
  df <- data.frame(x, y)
  df <- df[complete.cases(df), ]
  irr::icc(df, model = "twoway", type = "agreement", unit = "single")$value
}

text_base_size   <- 20
title_size       <- 20
axis_title_size  <- 20
axis_text_size   <- 18
annot_size       <- 18

lw_identity <- 0.8
lw_smooth   <- 1.0
lw_h0       <- 0.6
lw_bias     <- 1.0
lw_loa      <- 0.8

scatter_xlab <- "OEF (SJTU)"
scatter_ylab <- "OEF (ZJU)"
ba_xlab <- "Mean OEF"
ba_ylab <- "Difference"

common_ylim_sc <- c(0.25, 0.45)
common_xlim_sc <- c(0.25, 0.45)
common_ylim_ba <- c(-0.1, 0.1)
common_xlim_ba <- c(0.25, 0.45)

make_scatter <- function(df, x, y, color, title, filename,
                         xlim_range = common_xlim_sc,
                         ylim_range = common_ylim_sc,
                         xlab = scatter_xlab,
                         ylab = scatter_ylab){
  ct <- cor.test(df[[x]], df[[y]])
  r  <- unname(ct$estimate)
  p  <- ct$p.value
  icc_val <- calc_icc(df[[x]], df[[y]])
  x_span <- diff(xlim_range)
  y_span <- diff(ylim_range)
  
  p1 <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_abline(slope = 1, intercept = 0, color = "grey50",
                linetype = "dashed", linewidth = lw_identity) +
    geom_point(size = 3, color = color, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = color, fill = color,
                alpha = 0.18, linewidth = lw_smooth) +
    coord_equal() +
    coord_cartesian(xlim = xlim_range, ylim = ylim_range) +
    labs(title = title, x = xlab, y = ylab) +
    annotate("text",
             x = xlim_range[1] + 0.05 * x_span,
             y = ylim_range[2] - 0.1 * y_span,
             hjust = 0,
             size = annot_size / ggplot2::.pt,
             label = paste0("R=", round(r, 2),
                            "\nP=", signif(p, 2),
                            "\nICC=", round(icc_val, 2))) +
    theme_classic(base_size = text_base_size) +
    theme(
      plot.title = element_text(size = title_size, hjust = 0.5),
      axis.title = element_text(size = axis_title_size),
      axis.text  = element_text(color = "black", size = axis_text_size)
    )
  
  ggsave(file.path(outdir, filename), p1, width = 5.5, height = 5, dpi = 400)
}

make_ba <- function(df, x, y, color, title, filename,
                    xlim_range = common_xlim_ba,
                    ylim_range = common_ylim_ba,
                    xlab = ba_xlab,
                    ylab = ba_ylab){
  df2 <- df %>%
    transmute(
      mean = (.data[[x]] + .data[[y]]) / 2,
      diff = .data[[y]] - .data[[x]]
    ) %>%
    filter(complete.cases(.))
  
  bias  <- mean(df2$diff, na.rm = TRUE)
  sd_d  <- sd(df2$diff, na.rm = TRUE)
  loa_u <- bias + 1.96 * sd_d
  loa_l <- bias - 1.96 * sd_d
  
  p2 <- ggplot(df2, aes(x = mean, y = diff)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = lw_h0) +
    geom_hline(yintercept = bias, color = color, linewidth = lw_bias) +
    geom_hline(yintercept = loa_u, linetype = "dashed",
               color = "grey40", linewidth = lw_loa) +
    geom_hline(yintercept = loa_l, linetype = "dashed",
               color = "grey40", linewidth = lw_loa) +
    geom_point(size = 3, color = color, alpha = 0.9) +
    coord_cartesian(xlim = xlim_range, ylim = ylim_range) +
    labs(title = title, x = xlab, y = ylab) +
    theme_classic(base_size = text_base_size) +
    theme(
      plot.title = element_text(size = title_size, hjust = 0.5),
      axis.title = element_text(size = axis_title_size),
      axis.text  = element_text(color = "black", size = axis_text_size)
    )
  
  ggsave(file.path(outdir, filename), p2, width = 5.5, height = 5, dpi = 400)
}

site_rep1 <- dat %>%
  select(name, OEF_SJTU, OEF_ZJU)

make_scatter(site_rep1,
             "OEF_SJTU", "OEF_ZJU",
             "mediumpurple3",
             "Correlation between Sites",
             "Site_Correlation.png")

make_ba(site_rep1,
        "OEF_SJTU", "OEF_ZJU",
        "mediumpurple3",
        "Inter-site Bland–Altman",
        "Inter_site_BA.png")

site_rep1_long <- site_rep1 %>%
  pivot_longer(cols = c(OEF_SJTU, OEF_ZJU),
               names_to = "Site",
               values_to = "OEF") %>%
  mutate(Site = factor(Site,
                       levels = c("OEF_SJTU", "OEF_ZJU"),
                       labels = c("SJTU", "ZJU")))

p_line <- ggplot(site_rep1_long, aes(x = Site, y = OEF, group = name)) +
  geom_line(color = "grey30", linewidth = 0.9, alpha = 0.9) +
  geom_point(aes(color = name), size = 3, alpha = 0.95) +
  coord_cartesian(ylim = common_ylim_sc) +
  labs(title = "Inter-site Paired Plot", x = "", y = "OEF") +
  theme_classic(base_size = text_base_size) +
  theme(
    plot.title = element_text(size = title_size, hjust = 0.5),
    axis.title = element_text(size = axis_title_size),
    axis.text  = element_text(color = "black", size = axis_text_size),
    legend.position = "none"
  ) +
  scale_color_hue()

ggsave(file.path(outdir, "Inter_site_Paired_Line.png"),
       p_line, width = 5.5, height = 5, dpi = 400)
