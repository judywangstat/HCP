## Table 1 results
rm(list = ls())
if (exists("cl")) {
  stopCluster(cl)
  rm(cl)
}
gc()

library(ggplot2)
library(gridExtra)
library(rstudioapi)
main_path <- getSourceEditorContext()$path
main_dir <- dirname(main_path)


S_values <- c(1, 3, 5, 9, 11, 13, 15, 17, 19, 21)
n_values <- c(100, 300, 500)

results_matrix <- matrix(NA, nrow = length(S_values), ncol = length(n_values) * 3)
rownames(results_matrix) <- paste0("S=", S_values)
colnames(results_matrix) <- c(
  paste0("n=", n_values[1], "_coverage"), paste0("n=", n_values[1], "_length"), paste0("n=", n_values[1], "_sd"),
  paste0("n=", n_values[2], "_coverage"), paste0("n=", n_values[2], "_length"), paste0("n=", n_values[2], "_sd"),
  paste0("n=", n_values[3], "_coverage"), paste0("n=", n_values[3], "_length"), paste0("n=", n_values[3], "_sd")
)

for (i in seq_along(S_values)) {
  S <- S_values[i]
  for (j in seq_along(n_values)) {
    n <- n_values[j]
    
    assign("S_splitting", S, envir = .GlobalEnv)
    assign("n", n, envir = .GlobalEnv)
    
    hcp_out <- capture.output(source(file.path(main_dir, "HCP_marginal.R"), local = TRUE, chdir = TRUE))
    hcp_vals <-  as.numeric(strsplit(trimws(hcp_out), "\\s+")[[1]])
    
    if (length(hcp_vals) == 3) {
      results_matrix[i, (j - 1) * 3 + 1] <- hcp_vals[1]  # coverage
      results_matrix[i, (j - 1) * 3 + 2] <- hcp_vals[2]  # length
      results_matrix[i, (j - 1) * 3 + 3] <- hcp_vals[3]  # sd
    } else {
      warning(paste("Unexpected output format for S =", S, "and n =", n))
    }
  }
}



## plot coverage and length
cover_plot <- function(data, sample_size, y_var) {
  ggplot(data, aes(x = S)) + 
    geom_line(aes_string(y = y_var), color = 'blue') + 
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
    labs(x = "S", y = "Marginal Coverage", title = paste("n =", sample_size)) +
    ylim(0.8, 1.0) +
    theme_minimal(base_size = 7.5)+
    theme(plot.title = element_text(hjust = 0.5))
}

leng_plot <- function(data, sample_size, y_col, sd_col) {
  ggplot(data) +
    geom_line(aes_string(x = "S", y = y_col), color = "blue") +
    geom_ribbon(aes_string(x = "S", ymin = paste0(y_col, " - 1.645 * ", sd_col), ymax = paste0(y_col, " + 1.645 * ", sd_col)), 
                fill = "blue", alpha = 0.2) +
    labs(x = "S", y = "Marginal Length", title = paste("n =", sample_size)) +
    theme_minimal(base_size = 7.5) +
    theme(plot.title = element_text(hjust = 0.5))+
    ylim(13, 22)
}


results_df <- data.frame(
  S = rep(S_values, each = 1),
  n100_coverage = results_matrix[, 1],
  n100_length = results_matrix[, 2],
  n100_sd = results_matrix[, 3],
  n300_coverage = results_matrix[, 4],
  n300_length = results_matrix[, 5],
  n300_sd = results_matrix[, 6],
  n500_coverage = results_matrix[, 7],
  n500_length = results_matrix[, 8],
  n500_sd = results_matrix[, 9]
)

p1 <- cover_plot(results_df, 100, "n100_coverage")
p2 <- cover_plot(results_df, 300, "n300_coverage")
p3 <- cover_plot(results_df, 500, "n500_coverage")
p4 <- leng_plot(results_df, 100, "n100_length", "n100_sd")
p5 <- leng_plot(results_df, 300, "n300_length", "n300_sd")
p6 <- leng_plot(results_df, 500, "n500_length", "n500_sd")



plot_grid = grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)
ggsave(filename = file.path(main_dir, "sens_S_reproduce.pdf"), plot = plot_grid, width = 7, height = 5)
write.csv(results_matrix, file = file.path(main_dir, "final_results_reproduce.csv"), row.names = FALSE)
