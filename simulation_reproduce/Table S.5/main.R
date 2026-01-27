## Table S.4 results_all
rm(list = ls())
if (exists("cl")) {
  stopCluster(cl)
  rm(cl)
}
gc()

library(rstudioapi)
main_path <- getSourceEditorContext()$path
main_dir <- dirname(main_path)

sample_sizes <- c(100, 300, 500)
scenarios <- c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix")

n_rows <- length(sample_sizes) * length(scenarios)
results_all <- matrix(NA, nrow = n_rows, ncol = 4)
colnames(results_all) <- c("density_cov", "density_len", "residual_cov", "residual_len")

results_info <- data.frame(sample_size = rep(sample_sizes, each = length(scenarios)),
                           scenario = rep(scenarios, times = length(sample_sizes)))

row_index <- 1
for(n_val in sample_sizes) {
  for(scen in scenarios) {
    assign("n", n_val, envir = .GlobalEnv)
    assign("scenario", scen, envir = .GlobalEnv)
    
    out <- capture.output(source(file.path(main_dir,"HCP_compare_score.R"), local = TRUE, chdir = TRUE))
    vals <- as.numeric(strsplit(trimws(out), "\\s+")[[1]])
    
    
    results_all[row_index, ] <- vals
    row_index <- row_index + 1
  }
}

final_results <- cbind(results_info, results_all)
print(final_results)


write.csv(final_results, file = file.path(main_dir, "final_results_reproduce.csv"), row.names = FALSE)


