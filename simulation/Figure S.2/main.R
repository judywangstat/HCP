rm(list = ls())
if (exists("cl")) {
  stopCluster(cl)
  rm(cl)
}
gc()

library(ggplot2)
library(gridExtra)
library(rstudioapi)


# Get the directory of the current script
main_path <- getSourceEditorContext()$path
main_dir <- dirname(main_path)

# Function to run HCP_condition.R with different values of nclusters
run_HCP_condition <- function(nclusters) {
  assign("nclusters", nclusters, envir = .GlobalEnv)
  source(file.path(main_dir, "HCP_conditional.R"), local = FALSE, chdir = TRUE)
  return(output_matrix)
}

output_matrix_1 <- run_HCP_condition(nclusters = 1)
output_matrix_5 <- run_HCP_condition(nclusters = 5)


combined_matrix <- cbind(output_matrix_1, output_matrix_5)
colnames(combined_matrix) <- c("Mean_Cover_1", "SD_Cover_1", "Mean_Leng_1", "SD_Leng_1",
                               "Mean_Cover_5", "SD_Cover_5", "Mean_Leng_5", "SD_Leng_5")

print(combined_matrix)

# Generate index values for x-axis
index <- seq(-3, 3, length.out = nrow(combined_matrix))

# Create a data frame for plotting
df <- data.frame(
  index = index,
  Mean_Cover_1 = combined_matrix[, "Mean_Cover_1"],
  SD_Cover_1 = combined_matrix[, "SD_Cover_1"],
  Mean_Leng_1 = combined_matrix[, "Mean_Leng_1"],
  SD_Leng_1 = combined_matrix[, "SD_Leng_1"],
  Mean_Cover_5 = combined_matrix[, "Mean_Cover_5"],
  SD_Cover_5 = combined_matrix[, "SD_Cover_5"],
  Mean_Leng_5 = combined_matrix[, "Mean_Leng_5"],
  SD_Leng_5 = combined_matrix[, "SD_Leng_5"]
)


# Conditional Coverage 
cover_plot <- function(df, y_var, y_sd_var, title) {
  ggplot(df, aes(x = index)) + 
    geom_line(aes_string(y = y_var), color = 'blue') + 
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
    geom_ribbon(aes_string(ymin = paste(y_var, "- 1.645 *", y_sd_var), 
                           ymax = paste(y_var, "+ 1.645 *", y_sd_var)), 
                alpha = 0.2, fill = 'blue') + 
    labs(title = title, x = expression(x[2]), y = "Conditional Coverage") +
    theme_minimal(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    scale_y_continuous(limits = c(0.6, 1.08), breaks = seq(0.6, 1.0, by = 0.1)) # Set y-axis range and ticks
}

# Conditional Length 
leng_plot <- function(df, y_var, y_sd_var, title) {
  ggplot(df, aes(x = index)) + 
    geom_line(aes_string(y = y_var), color = 'blue') + 
    geom_ribbon(aes_string(ymin = paste(y_var, "- 3 * 1.645 *", y_sd_var), 
                           ymax = paste(y_var, "+ 3 * 1.645 *", y_sd_var)), 
                alpha = 0.2, fill = 'blue') + 
    labs(title = title, x = expression(x[2]), y = "Conditional Length") +
    ylim(0, 35) +
    theme_minimal(base_size = 8) +
    theme(plot.title = element_text(hjust = 0.5)) # Center the title
}

# Generate plots for both K=1 and K=5 cases
p1 <- cover_plot(df, "Mean_Cover_1", "SD_Cover_1", "K=1")
p2 <- leng_plot(df, "Mean_Leng_1", "SD_Leng_1", "K=1")
p3 <- cover_plot(df, "Mean_Cover_5", "SD_Cover_5", "K=5")
p4 <- leng_plot(df, "Mean_Leng_5", "SD_Leng_5", "K=5")
plot_grid = grid.arrange(p1, p3, p2, p4, ncol = 2)

# Save PDF plot and CSV file
ggsave(filename = file.path(main_dir, "condi_100_reproduce.pdf"), plot = plot_grid, width = 5.5, height = 4)
write.csv(combined_matrix, file = file.path(main_dir, "final_results_reproduce.csv"), row.names = FALSE)

