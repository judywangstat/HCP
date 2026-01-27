## Table 2 results
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

results_all <- matrix(NA, nrow = 6, ncol = 10)
colnames(results_all) <- c("Homo_cov", "Homo_len", 
                       "Heter_cov", "Heter_len", 
                       "Asym_cov", "Asym_len", 
                       "Bimo_cov", "Bimo_len", 
                       "Bimo-fix_cov", "Bimo-fix_len")

missing_types <- rep(c("20%", "50%"), each = length(sample_sizes))
results_info <- data.frame(missing = missing_types, sample_size = rep(sample_sizes, times = 2))

row_index <- 1
for(n_val in sample_sizes) {
  scenario_results <- numeric(10)  
  col_index <- 1
  for(scen in scenarios) {
    assign("n", n_val, envir = .GlobalEnv)
    assign("scenario", scen, envir = .GlobalEnv)
    
    out <- capture.output(source(file.path(main_dir, "HCP_simul20.R"), local = TRUE, chdir = TRUE))
    vals <- as.numeric(strsplit(trimws(out), "\\s+")[[1]])
    
    scenario_results[col_index]   <- vals[1]
    scenario_results[col_index+1] <- vals[2]
    col_index <- col_index + 2
  }
  results_all[row_index, ] <- scenario_results
  row_index <- row_index + 1
}

for(n_val in sample_sizes) {
  scenario_results <- numeric(10)
  col_index <- 1
  for(scen in scenarios) {
    assign("n", n_val, envir = .GlobalEnv)
    assign("scenario", scen, envir = .GlobalEnv)
    
    out <- capture.output(source(file.path(main_dir, "HCP_simul50.R"), local = TRUE, chdir = TRUE))
    vals <- as.numeric(strsplit(trimws(out), "\\s+")[[1]])
    
    scenario_results[col_index]   <- vals[1]
    scenario_results[col_index+1] <- vals[2]
    col_index <- col_index + 2
  }
  results_all[row_index, ] <- scenario_results
  row_index <- row_index + 1
}

final_results <- cbind(results_info, results_all)
print(final_results)
write.csv(final_results, file = file.path(main_dir, "final_results_reproduce.csv"), row.names = FALSE)
