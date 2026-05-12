
library(readxl)
library(dplyr)
library(ggplot2)
library(irr)
library(lme4)
library(lmerTest)

file_path <- "data/trust_caffeine_followup.xlsx"
outdir    <- "outputs/reproducibility_caffeine"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

base_font_size  <- 20
axis_title_size <- 20
axis_text_size  <- 18
title_size      <- 20
legend_text_size <- 18

line_width  <- 1.1
point_size  <- 2.2
point_alpha <- 0.85

reg_line_width <- 1.0
ci_alpha_fill  <- 0.35
identity_lty   <- "dashed"

w_curve   <- 5.5
h_curve   <- 5
w_scatter <- 5.5
h_scatter <- 5
dpi_out   <- 300

n_boot <- 2000
seed_boot <- 1234

xlim_agree <- c(0.15, 0.60)
ylim_agree <- c(0.15, 0.60)

dat <- read_excel(file_path) %>%
  as.data.frame() %>%
  mutate(Name = as.factor(Name))

stopifnot(all(c("Name", "Time", "OEF_Day1", "OEF_Day2") %in% names(dat)))

dat$Time_f <- factor(dat$Time)

m_day1 <- lmer(OEF_Day1 ~ Time_f + (1 | Name), data = dat)
p_time_day1 <- anova(m_day1)["Time_f", "Pr(>F)"]

m_day2 <- lmer(OEF_Day2 ~ Time_f + (1 | Name), data = dat)
p_time_day2 <- anova(m_day2)["Time_f", "Pr(>F)"]

subjects <- levels(dat$Name)
n_sub <- length(subjects)

pal12 <- c(
  "lightcoral",
  "skyblue2",
  "mediumpurple3",
  "#00A087",
  "salmon2",
  "#3C5488",
  "slategray4",
  "palegreen2",
  "orchid3",
  "mistyrose3"
)

palette_sub <- pal12[seq_len(n_sub)]
names(palette_sub) <- subjects

theme_pub <- theme_classic(base_size = base_font_size, base_family = "Arial") +
  theme(
    plot.title = element_text(size = title_size, face = "plain", hjust = 0.5),
    axis.title = element_text(size = axis_title_size),
    axis.text  = element_text(size = axis_text_size),
    legend.position = "none"
  )

fmt_p <- function(p){
  if(is.na(p)) return("NA")
  if(p < 0.001) return("< 0.001")
  sprintf("= %.3f", p)
}

orthogonal_fit <- function(x, y){
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if(length(x) < 2){
    return(list(intercept = NA_real_, slope = NA_real_))
  }
  
  xm <- mean(x)
  ym <- mean(y)
  
  Xc <- cbind(x - xm, y - ym)
  S  <- cov(Xc)
  
  eig <- eigen(S)
  v1  <- eig$vectors[, 1]
  
  if(abs(v1[1]) < .Machine$double.eps){
    return(list(intercept = NA_real_, slope = NA_real_))
  }
  
  slope <- v1[2] / v1[1]
  intercept <- ym - slope * xm
  
  list(intercept = intercept, slope = slope)
}

predict_orthogonal <- function(x, fit){
  fit$intercept + fit$slope * x
}

bootstrap_orthogonal_ci <- function(data, xvar, yvar, x_grid, n_boot = 2000, seed = 1234){
  set.seed(seed)
  
  dsub <- data %>%
    dplyr::select(all_of(c(xvar, yvar))) %>%
    dplyr::filter(is.finite(.data[[xvar]]), is.finite(.data[[yvar]]))
  
  n <- nrow(dsub)
  if(n < 5){
    return(data.frame(
      x = x_grid,
      y = NA_real_,
      y_low = NA_real_,
      y_high = NA_real_
    ))
  }
  
  fit0 <- orthogonal_fit(dsub[[xvar]], dsub[[yvar]])
  y0 <- predict_orthogonal(x_grid, fit0)
  
  boot_mat <- matrix(NA_real_, nrow = length(x_grid), ncol = n_boot)
  
  for(b in seq_len(n_boot)){
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    xb <- dsub[[xvar]][idx]
    yb <- dsub[[yvar]][idx]
    
    fitb <- orthogonal_fit(xb, yb)
    
    if(is.finite(fitb$slope) && is.finite(fitb$intercept)){
      boot_mat[, b] <- predict_orthogonal(x_grid, fitb)
    }
  }
  
  y_low  <- apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
  y_high <- apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  data.frame(
    x = x_grid,
    y = y0,
    y_low = y_low,
    y_high = y_high
  )
}

p1 <- ggplot(dat, aes(x = Time, y = OEF_Day1, group = Name, color = Name)) +
  geom_line(linewidth = line_width, alpha = 1) +
  geom_point(size = point_size, alpha = point_alpha) +
  scale_color_manual(values = palette_sub) +
  labs(
    x = "Time (min)",
    y = "OEF (Day 1)",
    title = "OEF after caffeine challenge (Day 1)"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = paste0("P ", fmt_p(p_time_day1)),
    hjust = -0.1, vjust = 1.2,
    family = "Arial",
    size = (base_font_size - 2) / ggplot2::.pt
  ) +
  coord_cartesian(ylim = c(0.15, 0.60)) +
  theme_pub

p2 <- ggplot(dat, aes(x = Time, y = OEF_Day2, group = Name, color = Name)) +
  geom_line(linewidth = line_width, alpha = 1) +
  geom_point(size = point_size, alpha = point_alpha) +
  scale_color_manual(values = palette_sub) +
  labs(
    x = "Time (min)",
    y = "OEF (Day 2)",
    title = "OEF after caffeine challenge (Day 2)"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = paste0("P ", fmt_p(p_time_day2)),
    hjust = -0.1, vjust = 1.2,
    family = "Arial",
    size = (base_font_size - 2) / ggplot2::.pt
  ) +
  coord_cartesian(ylim = c(0.15, 0.60)) +
  theme_pub

dat_icc <- dat %>%
  mutate(Target = paste(Name, Time, sep = "__")) %>%
  select(Target, OEF_Day1, OEF_Day2) %>%
  arrange(Target)

ct <- cor.test(dat$OEF_Day1, dat$OEF_Day2, method = "pearson")
R_val <- unname(ct$estimate)
p_val <- ct$p.value

icc_res <- irr::icc(
  dat_icc[, c("OEF_Day1", "OEF_Day2")],
  model = "twoway",
  type  = "agreement",
  unit  = "single"
)
ICC_val <- icc_res$value

x_grid <- seq(xlim_agree[1], xlim_agree[2], length.out = 200)
ci_df <- bootstrap_orthogonal_ci(
  data = dat,
  xvar = "OEF_Day1",
  yvar = "OEF_Day2",
  x_grid = x_grid,
  n_boot = n_boot,
  seed = seed_boot
)

ortho_fit0 <- orthogonal_fit(dat$OEF_Day1, dat$OEF_Day2)
ortho_slope <- ortho_fit0$slope
ortho_intercept <- ortho_fit0$intercept

label_text <- paste0(
  "R = ", sprintf("%.2f", R_val),
  "\nP ", fmt_p(p_val),
  "\nICC = ", sprintf("%.2f", ICC_val)
)

p3 <- ggplot(dat, aes(x = OEF_Day1, y = OEF_Day2, color = Name)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = palette_sub) +
  geom_abline(slope = 1, intercept = 0, linetype = identity_lty, linewidth = 0.8) +
  geom_ribbon(
    data = ci_df,
    aes(x = x, ymin = y_low, ymax = y_high),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = ci_alpha_fill
  ) +
  geom_line(
    data = ci_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color = "black",
    linewidth = reg_line_width
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = label_text,
    hjust = -0.2, vjust = 1.1,
    family = "Arial",
    size = (base_font_size - 2) / ggplot2::.pt
  ) +
  coord_cartesian(
    xlim = xlim_agree,
    ylim = ylim_agree
  ) +
  labs(
    x = "OEF (Day 1)",
    y = "OEF (Day 2)",
    title = "Agreement between Day 1 and 2"
  ) +
  theme_pub

ggsave(file.path(outdir, "Fig1_OEF_Day1_vs_Time.png"),
       p1, width = w_curve, height = h_curve, dpi = dpi_out)

ggsave(file.path(outdir, "Fig2_OEF_Day2_vs_Time.png"),
       p2, width = w_curve, height = h_curve, dpi = dpi_out)

ggsave(file.path(outdir, "Fig3_Agreement_Day1_vs_Day2_OrthogonalCI.png"),
       p3, width = w_scatter, height = h_scatter, dpi = dpi_out)

cat("Day 1 time effect p:", p_time_day1, "\n")
cat("Day 2 time effect p:", p_time_day2, "\n")
cat("Pearson R:", R_val, "\n")
cat("Pearson p:", p_val, "\n")
cat("ICC(2,1):", ICC_val, "\n")
cat("Orthogonal slope:", ortho_slope, "\n")
cat("Orthogonal intercept:", ortho_intercept, "\n")

p1
p2
p3
