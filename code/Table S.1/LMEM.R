
# Load necessary libraries
library(MASS)
library(stats)
library(doParallel)
library(doRNG)
library(lme4)
library(merTools)
source("~/Desktop/MAR_R_code/my_functions.R")


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

clusterEvalQ(cl, {
  library(MASS)
  library(stats)
  library(lme4)
  library(merTools)
  source("~/Desktop/MAR_R_code/my_functions.R")
})


results <- foreach(
  i = 1:num_simulations,
  .combine = rbind,
  .options.RNG = 123
) %dorng% {
  # results = matrix(0, nrow = num_simulations, ncol = 4)
  #   for(i in 1:num_simulations) {
  #   print(i)
  #   set.seed(i+100)
  
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
  
  # point_pred = predict(random_effect_model,
  #                      newdata = data_test,
  #                      allow.new.levels = TRUE)
  #
  # residual_var <- sigma(random_effect_model)^2
  # beta_cov_matrix <- vcov(random_effect_model)
  # random_effects_cov <- as.matrix(VarCorr(random_effect_model)$group)
  # fixed_effect_var = diag(X_test[[1]] %*% beta_cov_matrix %*% t(X_test[[1]]))
  # random_effect_var  =  diag(X_test[[1]] %*% random_effects_cov %*% t(X_test[[1]]))
  #
  # total_variance = residual_var + fixed_effect_var + random_effect_var
  #
  #
  # residual_sd <- sqrt(total_variance)
  # #residual_sd <- sigma(random_effect_model)
  # z_value <- qnorm(1 - alpha / 2)
  # lower_bound <- point_pred - z_value * residual_sd
  # upper_bound <- point_pred + z_value * residual_sd
  #
  # pred_inter <- data.frame(
  #   Y_true = Y_test[[1]],
  #   Prediction = point_pred,
  #   Lower_PI_90 = lower_bound,
  #   Upper_PI_90 = upper_bound
  # )
  
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
