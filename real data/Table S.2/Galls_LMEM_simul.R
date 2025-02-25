
# Load necessary libraries
library(MASS)
library(stats)
library(grf)
library(doParallel)
library(doRNG)
library(quantreg)
library(grf)
library(lme4)
library(merTools)
library(randomForest)

library(rstudioapi)
script_path <- getSourceEditorContext()$path
script_dir <- dirname(script_path)
parent_dir <- dirname(script_dir)
source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)


### General parameters ###
alpha = 0.1/4
n_grid = 200
n_test <- 1

### Importing data ###
data <- read.delim(file.path(script_dir,"gallstones.txt"))
ALL <- data
Time <- ALL$time
Treat <- ALL$trt

T1 <- ifelse(Time == 12, 1, 0)
T2 <- ifelse(Time == 20, 1, 0)
T3 <- ifelse(Time == 24, 1, 0)
x <- cbind(T1, T2, T3, Treat)
y <- as.numeric(ALL$y)
id <- ALL$id

n <- length(unique(id))
X_list <- vector("list", n)
Y_list <- vector("list", n)
Delta_list <- vector("list", n)

for (ii in 1:n) {
  X_list[[ii]] <- x[id == ii,]
  Y_list[[ii]] <- y[id == ii]
  Delta_list[[ii]] <- ifelse(is.na(y[id == ii]), 0, 1)  
}

### Imputing missing values ###
set.seed(11111)
X_matrix <- do.call(rbind, X_list)
Y_vector <- unlist(Y_list)

non_na_index <- !is.na(Y_vector)
na_index <- is.na(Y_vector)

m_rf_model <- randomForest(X_matrix[non_na_index, ], Y_vector[non_na_index])
m_X_pred <- predict(m_rf_model, X_matrix[na_index, ])


quantile_forest_r <- quantile_forest(X_matrix[non_na_index, ], Y_vector[non_na_index], quantiles = c(0.25, 0.75))
quantiles_pred <- predict(quantile_forest_r, X_matrix[na_index, ])$predictions
r_X_pred <- quantiles_pred[, 2] - quantiles_pred[, 1]
e <- rnorm(length(m_X_pred))
Y_pred <- m_X_pred +  1.2* r_X_pred * e
Y_vector[na_index] <- Y_pred
Y_list_imputed <- split(Y_vector, rep(1:length(X_list), each = 4))



num_cores <- detectCores()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterExport(cl, "parent_dir")

clusterEvalQ(cl, {
  library(MASS)
  library(stats)
  library(grf)
  library(lme4)
  library(merTools)
  library(quantreg)
  source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)
})

results <- foreach(
  test_id = seq_len(n),  
  .combine = rbind,  
  .options.RNG = 123  
) %dorng% {
  
  
  ### Leave one out test data ###
  X_test <- X_list[[test_id]]
  Y_test <- Y_list_imputed[[test_id]]
  Delta_test <- Delta_list[[test_id]]
  n_simultaneous = length(Y_test)
  
  X_sample <- X_list[-test_id]
  Y_sample <- Y_list[-test_id]
  Delta_sample <- Delta_list[-test_id]
  
  n_sample = length(Y_sample)
  group <- factor(rep(1:n_sample, each = 4))
  
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
  random_effect_model <- lmer(Y ~ 1 + X1 + X2 + X3 + X4 + ( 1 | group), data = data)
  
  data_test <-data.frame(group = rep(0, 4))
  data_test[X_columns] <- X_test
  
  pred_intervals <- predictInterval(random_effect_model, newdata = data_test, level = 1-alpha, n.sims = 1000, which = "full")
  
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







