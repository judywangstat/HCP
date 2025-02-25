## Table S.1 results
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
n_comb <- length(sample_sizes) * length(scenarios)
results_all <- matrix(NA, nrow = n_comb, ncol = 8)
colnames(results_all) <- c("HCP_cov", "HCP_len", "DWR_cov", "DWR_len",
                       "LC_cov", "LC_len", "LMEM_cov", "LMEM_len")


results_info <- data.frame(sample_size = rep(sample_sizes, each = length(scenarios)),
                           scenario = rep(scenarios, times = length(sample_sizes)))

row_index <- 1
for(n_val in sample_sizes) {
  for(scen in scenarios) {
    assign("n", n_val, envir = .GlobalEnv)
    assign("scenario", scen, envir = .GlobalEnv)
  
    hcp_out <- capture.output(source(file.path(main_dir, "HCP_marginal.R"), local = TRUE, chdir = TRUE))
    hcp_vals <-  as.numeric(strsplit(trimws(hcp_out), "\\s+")[[1]])
    
    dwr_out <- capture.output(source(file.path(main_dir, "DWR.R"), local = TRUE, chdir = TRUE))
    dwr_vals <- as.numeric(strsplit(trimws(dwr_out), "\\s+")[[1]])
    
    lc_out <- capture.output(source(file.path(main_dir, "LC.R"), local = TRUE, chdir = TRUE))
    lc_vals <- as.numeric(strsplit(trimws(lc_out), "\\s+")[[1]])
    
    lmem_out <- capture.output(source(file.path(main_dir, "LMEM.R"), local = TRUE, chdir = TRUE))
    lmem_vals <- as.numeric(strsplit(trimws(lmem_out), "\\s+")[[1]])
    
    results_all[row_index, ] <- c(hcp_vals[1], hcp_vals[2],
                              dwr_vals[1], dwr_vals[2],
                              lc_vals[1], lc_vals[2],
                              lmem_vals[1], lmem_vals[2])
    
    row_index <- row_index + 1
  }
}

final_results <- cbind(results_info, results_all)
print(final_results)

write.csv(final_results, file = file.path(main_dir, "final_results_reproduce.csv"), row.names = FALSE)





