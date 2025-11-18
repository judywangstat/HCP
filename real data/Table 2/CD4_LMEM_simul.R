
# Load necessary libraries
library(MASS)
library(stats)
library(doParallel)
library(doRNG)
library(lme4)
library(merTools)

library(rstudioapi)
script_path <- getSourceEditorContext()$path
script_dir <- dirname(script_path)
parent_dir <- dirname(script_dir)
source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)


### General parameters ###
scenario = NULL
#missing <- "20"  #c("20", "50")

### Importing data ###
cd4 <- read.table(file.path(script_dir,"CD4_data.txt"), header=TRUE, quote="\"")
x_vars <- cd4[, c("time", "age", "smoke", "drug", "partners", "cesd")]
y <- cd4$CD4
id <- cd4$id
n <- length(unique(id))
unid <- unique(id)

Ti <- integer(n)
X_list <- vector("list", n)
Y_list <- vector("list", n)

for (ii in seq_len(n)) {
  idx <- id == unid[ii]
  Ti[ii] <- sum(idx)   
  X_list[[ii]] <- as.matrix(x_vars[idx, ])
  Y_list[[ii]] <- y[idx]
}


### Generate missing indicator ###
set.seed(123)
beta <- switch(missing,
               "20" = c(7, 0, 0, -3, 0, 0, 0),
               "50" = c(1, 0, 0, -3, 0, 0, 0),
               stop("Invalid missing value")
)
Delta_list <- vector("list", n)
for (ii in seq_len(n)) {
  X <- X_list[[ii]]
  X <- as.matrix(X)
  Delta_list[[ii]] <- generate_missing_indicator(X, beta)
}



num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterExport(cl, "parent_dir")

clusterEvalQ(cl, {
  library(MASS)
  library(stats)
  library(lme4)
  library(merTools)
  source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)
})

results <- foreach(
  test_id = seq_len(n),  
  .combine = rbind,  
  .options.RNG = 123  
) %dorng% {


  ### Leave one out test data ###
  Ti_test <- Ti[test_id]
  X_test <- X_list[[test_id]]
  Y_test <- Y_list[[test_id]]
  Delta_test <- Delta_list[[test_id]]
  n_simultaneous = length(Y_test)
  alpha = 0.1/n_simultaneous
  
  
  Ti_sample <- Ti[-test_id]
  X_sample <- X_list[-test_id]
  Y_sample <- Y_list[-test_id]
  Delta_sample <- Delta_list[-test_id]
  n_sample = length(Y_sample)
  group <- generate_group(Ti_sample)
  
  
  
  ## using observed data
  Y_whole_sample <- unlist(Y_sample)
  Delta_sample <- unlist(Delta_sample)
  X_whole_sample <- do.call(rbind, X_sample)
  Y_obs_training <- Y_whole_sample[Delta_sample == 1]
  X_obs_training <- X_whole_sample[Delta_sample == 1, ]
  group_obs <- group[Delta_sample == 1]
  
  
  num_features <- ncol(X_obs_training)
  X_columns <- paste0("X", 1:num_features)
  data <- data.frame(
    group = group_obs,
    Y = Y_obs_training
  )
  data[X_columns] <- X_obs_training
  random_effect_model <- lmer(Y ~ 1 + X1 + X2 + X3 + X4 + X5 + X6 + (1  | group), data = data)
  
  
  
  data_test <-data.frame(group = rep(0, Ti_test))
  data_test[X_columns] <- X_test
  
  pred_intervals <- predictInterval(random_effect_model, newdata = data_test, level = 1- alpha, n.sims = 1000, which = "full")
  
  pred_inter <- data.frame(
    Y_true = Y_test,
    Prediction = pred_intervals$fit,
    Lower_PI_90 = pred_intervals$lwr,
    Upper_PI_90 = pred_intervals$upr
  )
  
  
  return(pred_inter)
}

stopCluster(cl)
covered <- (results$Y_true >= results$Lower_PI_90) &
  (results$Y_true <= results$Upper_PI_90)
covered_groups <- split(covered, ceiling(seq_along(covered) / 4))
group_coverage <- sapply(covered_groups, all)
interval_length <- results$Upper_PI_90 - results$Lower_PI_90


res <- c(mean(group_coverage), mean(interval_length))
cat(paste(res, collapse = " "), "\n")




