rm(list = ls())
if (exists("cl")) {
  stopCluster(cl)
  rm(cl)
}
gc()

library(rstudioapi)
main_path <- getSourceEditorContext()$path
main_dir <- dirname(main_path)

missing_settings  <- c("20", "50")
prediction_types  <- c("pointwise", "simul")
methods           <- c("HCP", "DWR", "LC", "LMEM")


final_table <- matrix(NA, nrow = 8, ncol = length(methods))
colnames(final_table) <- methods
rownames(final_table) <- c("pointwise_20_coverage", "pointwise_20_length",
                           "pointwise_50_coverage", "pointwise_50_length",
                           "simul_20_coverage",     "simul_20_length",
                           "simul_50_coverage",     "simul_50_length")



row_counter <- 1
for(pred in prediction_types) {
  for(miss in missing_settings) {
    assign("missing", miss, envir = .GlobalEnv)
    

    for(i in seq_along(methods)) {
      method <- methods[i]
      if (pred == "pointwise") {
        script_file <- paste0("CD4_", method, ".R")
      } else {
        script_file <- paste0("CD4_", method, "_simul.R")
      }
      
      full_script_path <- file.path(main_dir, script_file)
      out <- capture.output(source(full_script_path, local = TRUE, chdir = TRUE))
      vals <- as.numeric(strsplit(trimws(out), "\\s+")[[1]])
      
      # row_counter:     coverage
      # row_counter + 1: length
      final_table[row_counter, i]   <- vals[1]
      final_table[row_counter + 1, i] <- vals[2]
    }
    row_counter <- row_counter + 2
  }
}


print(final_table)
write.csv(final_table, file = file.path(main_dir, "final_results_reproduce.csv"), row.names = TRUE)




