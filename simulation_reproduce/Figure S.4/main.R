rm(list = ls())
if (exists("cl")) {
  stopCluster(cl)
  rm(cl)
}
gc()

library(ggplot2)
library(viridis)
library(gridExtra)
library(rstudioapi)


# Get the directory of the current script
main_path <- getSourceEditorContext()$path
main_dir <- dirname(main_path)

# Function to run HCP_local.R with different values of nclusters
run_HCP_local <- function(nclusters) {
  assign("nclusters", nclusters, envir = .GlobalEnv)
  source(file.path(main_dir, "HCP_local.R"), local = FALSE, chdir = TRUE)
  return(final_matrix)
}

output_matrix_1 <- run_HCP_local(nclusters = 1)
output_matrix_10 <- run_HCP_local(nclusters = 10)
combined_matrix <- rbind(output_matrix_1, output_matrix_10)
rownames(combined_matrix) <- c("Coverage_K1", "Length_K1", "Coverage_K10", "Length_K10")

x2_breaks <- seq(-2.5, 2.5, length.out = 6)
x3_breaks <- seq(-2.5, 2.5, length.out = 6)
grid_data <- expand.grid(x2 = x2_breaks, x3 = x3_breaks)

grid_data$coverage_prob_K1 <- output_matrix_1[1, ]  # First row: Coverage for K=1
grid_data$interval_length_K1 <- output_matrix_1[2, ]    # Second row: Length for K=1
grid_data$coverage_prob_K10  <- output_matrix_10[1, ] # Third row: Coverage for K=10
grid_data$interval_length_K10 <- output_matrix_10[2, ]   # Fourth row: Length for K=10


limits_range_cover <- range(c(grid_data$coverage_prob_K1, grid_data$coverage_prob_K10), na.rm = TRUE)
limits_range_leng <- range(c(grid_data$interval_length_K1, grid_data$interval_length_K10), na.rm = TRUE)



p1 <- ggplot(grid_data, aes(x = x2, y = x3, fill = coverage_prob_K1)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "D", direction = 1, limits = limits_range_cover) +
  theme_minimal(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = expression(x[2]), y = expression(x[3]), title = "Local Coverage, K=1", fill = " ")

p2 <- ggplot(grid_data, aes(x = x2, y = x3, fill = interval_length_K1)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "D", direction = 1, limits = limits_range_leng) +
  theme_minimal(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = expression(x[2]), y = expression(x[3]), title = "Local Length, K=1", fill = "")

p3 <- ggplot(grid_data, aes(x = x2, y = x3, fill = coverage_prob_K10)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "D", direction = 1, limits = limits_range_cover) +
  theme_minimal(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = expression(x[2]), y = expression(x[3]), title = "Local Coverage, K=10", fill = " ")

p4 <- ggplot(grid_data, aes(x = x2, y = x3, fill = interval_length_K10)) +
  geom_tile(color = "white") +
  scale_fill_viridis(option = "D", direction = 1, limits = limits_range_leng) +
  theme_minimal(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = expression(x[2]), y = expression(x[3]), title = "Local Length, K=10", fill = "")


plot_grid = grid.arrange(p1, p3, p2, p4, nrow = 2, ncol = 2)

ggsave(filename = file.path(main_dir, "local_300_reproduce.pdf"), plot = plot_grid, width = 7, height = 5)
write.csv(combined_matrix, file = file.path(main_dir, "final_results_reproduce.csv"), row.names = FALSE)
