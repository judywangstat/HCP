## =========================================================
rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(ggh4x); library(rstudioapi)
})

## Detect script directory
script_path <- tryCatch(getSourceEditorContext()$path, error = function(e) "")
script_dir  <- if (nzchar(script_path)) dirname(script_path) else getwd()

## 1) Read results from five methods
DWR_final_mat    <- read.csv(file.path(script_dir, "DWR_final_mat.csv"),    stringsAsFactors = FALSE)
LC_final_mat     <- read.csv(file.path(script_dir, "LC_final_mat.csv"),     stringsAsFactors = FALSE)
LMEM_final_mat   <- read.csv(file.path(script_dir, "LMEM_final_mat.csv"),   stringsAsFactors = FALSE)
Oracle_final_mat <- read.csv(file.path(script_dir, "Oracle_final_mat.csv"), stringsAsFactors = FALSE)
HCP_final_mat    <- read.csv(file.path(script_dir, "HCP_final_mat.csv"),    stringsAsFactors = FALSE)

## 2) Settings and label definitions
settings_order  <- c("Homo","Heter","Asym","Bimo","Bimo.fix")
setting_labels  <- c(Homo="Homo", Heter="Heter", Asym="Asym", Bimo="Bimo", "Bimo.fix"="Bimo-fix")
method_order    <- c("HCP","DWR","LC","LMEM","Oracle")  # desired display order

## 3) Compute grouped means (K groups)
group_mean <- function(x, K = 50, seed = 123) {
  stopifnot(is.numeric(x))
  N <- length(x); stopifnot(N %% K == 0)
  set.seed(seed)
  gid <- sample(rep(1:K, each = N / K))
  as.numeric(tapply(x, gid, mean))
}

## 4) Construct long-format data for plotting
method_list <- list(
  "HCP"    = HCP_final_mat,
  "DWR"    = DWR_final_mat,
  "LC"     = LC_final_mat,
  "LMEM"   = LMEM_final_mat,
  "Oracle" = Oracle_final_mat
)

K <- 50
res_list <- list()

for (mtd in names(method_list)) {
  df <- method_list[[mtd]]
  need <- as.vector(rbind(paste0(settings_order, "_cov"),
                          paste0(settings_order, "_len")))
  stopifnot(all(need %in% names(df)))
  
  for (s in settings_order) {
    ## Coverage
    res_list[[length(res_list)+1]] <- data.frame(
      Setting = s, Metric = "Coverage", Method = mtd,
      Group = seq_len(K), Value = group_mean(df[[paste0(s,"_cov")]], K=K),
      stringsAsFactors = FALSE
    )
    ## Length
    res_list[[length(res_list)+1]] <- data.frame(
      Setting = s, Metric = "Length", Method = mtd,
      Group = seq_len(K), Value = group_mean(df[[paste0(s,"_len")]], K=K),
      stringsAsFactors = FALSE
    )
  }
}

## Final plotting dataframe
plot_df <- bind_rows(res_list) %>%
  mutate(
    Setting_raw = factor(Setting, levels = settings_order),
    Setting_lab = factor(Setting_raw, levels = settings_order,
                         labels = setting_labels[settings_order]),
    Metric      = factor(Metric, levels = c("Coverage","Length")),
    Method      = factor(Method, levels = method_order)
  )

## Vertical reference line for coverage = 0.9
vline_df <- plot_df %>%
  distinct(Setting_lab) %>%
  mutate(
    Metric = factor("Coverage", levels = levels(plot_df$Metric)),
    xint = 0.9
  )

desired_top_to_bottom <- c("HCP","DWR","LC","LMEM","Oracle")

## Main plot: horizontal boxplots + mean squares + reference lines
p_all <- ggplot(plot_df, aes(x = Value, y = Method)) +
  scale_y_discrete(limits = rev(desired_top_to_bottom)) +
  
  ## Boxplots
  geom_boxplot(width = 0.5, outlier.shape = NA,
               colour = "gray20", fill = "white", linewidth = 0.4) +
  
  ## Mean markers
  stat_summary(fun = mean, geom = "point",
               shape = 22, size = 1.2, stroke = 0.5, fill = "white") +
  
  ## Reference line only for Coverage metric
  geom_vline(data = vline_df, aes(xintercept = xint),
             colour = "red3", linewidth = 0.55) +
  
  ## Facet: 5 settings × 2 metrics, with independent x-scale
  ggh4x::facet_grid2(
    rows  = vars(Setting_lab),
    cols  = vars(Metric),
    scales = "free_x",
    independent = "x",
    switch = "x"        # column strip at bottom; row strip at right
  ) +
  
  ## Base theme
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    
    strip.placement     = "outside",
    strip.text.x        = element_text(face = "bold"),
    strip.text.y        = element_text(size = 10, face = "bold", angle = 270),
    axis.text           = element_text(size = 8),
    axis.text.y         = element_text(color = "black"),
    
    strip.background.y  = element_rect(fill = "grey88", colour = NA),
    strip.background.x  = element_blank(),
    
    axis.title.x        = element_blank(),
    axis.title.y        = element_text(size = 12),
    panel.spacing.x     = unit(1.2, "lines"),
    panel.spacing.y     = unit(1.2, "lines")
  ) +
  
  ## Scale for Coverage (Length uses free scale)
  ggh4x::facetted_pos_scales(
    x = list(
      Metric == "Coverage" ~ scale_x_continuous(
        limits = c(0.45, 1.0),
        breaks = seq(0.4, 1.0, by = 0.1),
        labels = function(z) sprintf("%.2f", z),
        expand = expansion(mult = c(0.00, 0.02))
      )
    )
  )

print(p_all)
