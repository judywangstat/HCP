
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


# General parameters
num_simulations <- 1000
#n <- 100  # Number of subjects/groups
n_test <- 1
m <- rep(5, n)  # Number of observations per subject/group
d <- 4  # Dimension of X and theta
#scenario = "Bimo-fix"  # c("Homo", "Heter", "Asym", "Bimo", "Bimo-fix")
alpha = 0.1
missing <- "50"  #c("20", "50")


# Parameters for distributions
theta_mu <- rep(2, d)
theta_sigma <- diag(1, d)
X_mu <- rep(0, d)
X_sigma <- diag(1, d)
if (missing == "20") {
  beta <- c(3, 0, 2, 2, 2)
} else if (missing == "50") {
  beta <- c(0, 0, 2, 2, 2)
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
  i = 1:num_simulations,
  .combine = rbind,
  .options.RNG = 123
) %dorng% {

  
  ###  generate data  ###
  data = generate_data(n, m, d, scenario, theta_mu, theta_sigma, X_mu, X_sigma, beta)
  ###  generate test data  ###
  test_index = single_split(n, n_test)
  X_test <- data$X[test_index$indices_part1]
  Y_test <- data$Y[test_index$indices_part1]
  Delta_test <- data$Delta[test_index$indices_part1]
  X_sample <- data$X[test_index$indices_part2]
  Y_sample <- data$Y[test_index$indices_part2]
  Delta_sample <- data$Delta[test_index$indices_part2]
  
  n_sample = length(Y_sample)
  group <- factor(rep(1:n_sample, each = m[1]))
  
  Y_whole_sample <- unlist(Y_sample)
  Delta_sample <- unlist(Delta_sample)
  X_whole_sample <- do.call(rbind, X_sample)
  
  ## using observed data
  Y_obs_training <- Y_whole_sample[Delta_sample == 1]
  X_obs_training <- X_whole_sample[Delta_sample == 1, ]
  group_obs <- group[Delta_sample == 1]
  
  num_features <- ncol(X_obs_training)
  X_columns <- paste0("X", 1:num_features)
  data <- data.frame(group = group_obs, Y = Y_obs_training)
  data[X_columns] <- X_obs_training
  
  
  random_effect_model <- lmer(Y ~ 1 + X1 + X2 + X3 + X4 + (0 + X3 + X4 |
                                                             group), data = data)
  
  data_test <- data.frame(group = rep(0, 5))
  data_test[X_columns] <- X_test[[1]]

  pred_intervals <- predictInterval(
    random_effect_model,
    newdata = data_test,
    level = 0.9,
    n.sims = 1000,
    which = "random"
  )
  
  pred_inter <- data.frame(
    Y_true = Y_test[[1]],
    Prediction = pred_intervals$fit,
    Lower_PI_90 = pred_intervals$lwr,
    Upper_PI_90 = pred_intervals$upr
  )
  
  return(pred_inter)
}


stopCluster(cl)


covered <- (results$Y_true >= results$Lower_PI_90) &
  (results$Y_true <= results$Upper_PI_90)
interval_length <- results$Upper_PI_90 - results$Lower_PI_90
res <- c(mean(covered), mean(interval_length))
cat(paste(res, collapse = " "), "\n")
