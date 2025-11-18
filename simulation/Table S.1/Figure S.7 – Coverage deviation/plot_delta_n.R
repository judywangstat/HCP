rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(rstudioapi)
  library(ggplot2)
  library(ggh4x)     # facet_grid2 for row-right & column-bottom style strips
  library(dplyr)
  library(tidyr)
})

## 0) Script directory
script_path <- tryCatch(getSourceEditorContext()$path, error = function(e) "")
script_dir  <- if (nzchar(script_path)) dirname(script_path) else getwd()

## 1) Sample sizes and file definitions
n_values <- c(100, 200, 300, 500, 800, 1200, 1500)

files <- c(
  Homo  = "Homo_results.csv",
  Heter = "Heter_results.csv",
  Asym  = "Asym_results.csv",
  Bimo  = "Bimo_results.csv",
  `Bimo-fix` = "Bimo-fix_results.csv"
)

## 2) Wilson confidence interval
wilson_ci <- function(k, n, conf.level = 0.95){
  if (n == 0) return(c(mean = NA, low = NA, high = NA))
  z <- qnorm(1 - (1 - conf.level)/2)
  p <- k / n
  denom  <- 1 + z^2/n
  center <- (p + z^2/(2*n)) / denom
  half   <- (z / denom) * sqrt(p*(1-p)/n + z^2/(4*n^2))
  c(mean = center, low = pmax(0, center - half), high = pmin(1, center + half))
}

## 3) Summarize a single result file into aggregated statistics
summarize_one <- function(dat){
  out <- data.frame(
    n = n_values,
    mean_cov = NA_real_,
    cov_low  = NA_real_,
    cov_high = NA_real_,
    mean_len = NA_real_,
    se_len   = NA_real_
  )
  for (i in seq_along(n_values)){
    n  <- n_values[i]
    cv <- dat[[paste0("cov_n", n)]]
    lv <- dat[[paste0("len_n", n)]]
    k     <- sum(cv == 1, na.rm = TRUE)
    n_eff <- sum(!is.na(cv))
    ci    <- wilson_ci(k, n_eff, 0.95)
    
    out$mean_cov[i] <- ci["mean"]
    out$cov_low[i]  <- ci["low"]
    out$cov_high[i] <- ci["high"]
    out$mean_len[i] <- mean(lv, na.rm = TRUE)
    out$se_len[i]   <- sd(lv,   na.rm = TRUE) / sqrt(sum(!is.na(lv)))
  }
  out$diff_cov <- abs(out$mean_cov - 0.9)
  out
}

## Unified plot theme
theme_mh <- function(base_size = 11){
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
      
      axis.text        = element_text(size = base_size - 3, color = "black"),
      axis.title.x     = element_text(size = base_size + 1, margin = margin(t = 8)),
      axis.title.y     = element_blank(),
      
      plot.title       = element_text(size = base_size + 3, face = "bold",
                                      hjust = 0, margin = margin(b = 8)),
      plot.title.position = "plot",
      plot.margin      = margin(8, 8, 8, 8),
      
      strip.placement    = "outside",
      strip.text.y       = element_text(size = 10, face = "bold", angle = 270),
      strip.background.y = element_rect(fill = "grey88", colour = NA),
      strip.text.x       = element_blank(),          # hide top strip labels
      strip.background.x = element_blank(),          # hide strip background
      strip.switch.pad.grid = unit(0, "mm"),
      
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(1.2, "lines")
    )
}

## 4) Read and merge all result files into long-format dataframe
setting_levels <- c("Homo","Heter","Asym","Bimo","Bimo-fix")

all_df <- bind_rows(lapply(names(files), function(st){
  fpath <- file.path(script_dir, files[[st]])
  dat   <- read.csv(fpath, stringsAsFactors = FALSE)
  df    <- summarize_one(dat)
  df$Setting <- st
  df
})) %>%
  mutate(
    Setting = factor(Setting, levels = setting_levels)
  ) %>%
  pivot_longer(
    cols = c(diff_cov, mean_len),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric,
                    levels = c("diff_cov","mean_len"),
                    labels = c("Coverage deviation", "Average length"))
  )

## 5) Baseline and shaded band for “Coverage deviation”
hline_df <- data.frame(
  Metric = factor("Coverage deviation", levels = levels(all_df$Metric)),
  yint = 0
)

band_df <- data.frame(
  Metric = factor("Coverage deviation", levels = levels(all_df$Metric)),
  xmin   = -Inf, xmax = Inf,
  ymin   = 0, ymax = 0.02
)

hline_df1 <- data.frame(
  Metric = factor("Coverage deviation", levels = levels(all_df$Metric)),
  yint = 0.02
)

## 6) Main plot (using ggh4x for row-right & column-bottom layout)
p_all <- ggplot(all_df, aes(x = n, y = Value, group = Setting)) +
  geom_rect(
    data = band_df,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "grey75", alpha = 0.4
  ) +
  geom_line(linewidth = 0.6, lineend = "round") +
  geom_point(shape = 21, size = 1.4, stroke = 0.5, fill = "white") +
  scale_x_continuous(breaks = c(500, 1000, 1500, 2000, 2500, 3000)) +
  labs(x = "Sample size (n)") +
  
  geom_hline(
    data = hline_df,
    aes(yintercept = yint),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_hline(
    data = hline_df1,
    aes(yintercept = yint),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  ggh4x::facet_grid2(
    rows  = vars(Setting),
    cols  = vars(Metric),
    scales = "free_y",
    independent = "y",
    switch = "x",
    labeller = labeller(
      Metric = as_labeller(
        c("Coverage deviation" = "bold(Coverage~deviation)",
          "Average length" = "bold(Average~length)"),
        default = label_parsed
      )
    )
  ) +
  ggh4x::facetted_pos_scales(
    y = list(
      Metric == "Coverage deviation" ~ scale_y_continuous(
        limits = c(-0.02, 0.10),
        breaks = seq(0, 0.10, by = 0.03),
        expand = expansion(mult = c(0.01, 0.03))
      ),
      Metric == "Average length" ~ scale_y_continuous(
        expand = expansion(mult = c(0.3, 0.5))
      )
    )
  ) +
  theme_mh()

print(p_all)


library(cowplot)

## Headers for two-column layout
hdr1 <- ggdraw(clip = "off") +
  draw_label(expression(bold(~~~~~Coverage~deviation)),
             y = 0.9, vjust = 1, size = 10) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

hdr2 <- ggdraw(clip = "off") +
  draw_label(expression(bold(Average~length)),
             y = 0.9, vjust = 1, size = 10) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

top_row <- plot_grid(hdr1, hdr2, ncol = 2)

plot_grid(
  top_row,
  p_all + theme(plot.margin = margin(t = 2, r = 8, b = 8, l = 8)),
  ncol = 1,
  rel_heights = c(0.025, 1)
)
