rm(list = ls())
if (exists("cl")) {
  stopCluster(cl)
  rm(cl)
}
gc()

library(rstudioapi)
main_path <- getSourceEditorContext()$path
main_dir <- dirname(main_path)

# Define prediction types and methods
prediction_types <- c("pointwise", "simul")
metrics <- c("coverage", "length")
methods <- c("HCP", "DWR", "LC", "LMEM")

# Initialize result table
final_table <- matrix(NA, nrow = 4, ncol = length(methods))
colnames(final_table) <- methods
rownames(final_table) <- c("pointwise_coverage", "pointwise_length", 
                           "simul_coverage", "simul_length")

# Row counter for the result table
row_counter <- 1
for (pred in prediction_types) {
  for (metric in metrics) {
    for (i in seq_along(methods)) {
      method <- methods[i]
      
      # Construct script filename
      if (pred == "pointwise") {
        script_file <- paste0("Galls_", method, ".R")
      } else {
        script_file <- paste0("Galls_", method, "_simul.R")
      }
      
      full_script_path <- file.path(main_dir, script_file)
      out <- capture.output(source(full_script_path, local = TRUE, chdir = TRUE))
      vals <- as.numeric(strsplit(trimws(out), "\\s+")[[1]])
      
      # Store values in the final table
      final_table[row_counter, i] <- vals[ifelse(metric == "coverage", 1, 2)]
    }
    row_counter <- row_counter + 1
  }
}

# Print and save results
print(final_table)
write.csv(final_table, file = file.path(main_dir, "final_results_reproduce.csv"), row.names = TRUE)



