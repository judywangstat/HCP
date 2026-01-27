
# Load necessary libraries
library(MASS)
library(stats)
library(grf)
library(doParallel)
library(doRNG)
library(quantreg)

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
#scenario = "Heter"  # "Homo", "Heter", "Asym", "Bimo", "Bimo-fix"
alpha = 0.1
n_grid = 200
quantile_model = "linear"   # "linear", "grf"
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



num_cores <- detectCores() - 3
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterExport(cl, "parent_dir")

clusterEvalQ(cl, {
  library(MASS)
  library(stats)
  library(grf)
  library(quantreg)
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
  
  lo_ygrid = min(unlist(Y_sample))
  up_ygrid = max(unlist(Y_sample))
  y_grid <- seq(from = lo_ygrid,
                to = up_ygrid,
                length.out = n_grid)
  
  
  ### Given test data
  Y_test = Y_test[[1]]
  X_test = X_test[[1]]
  n_simultaneous = length(Y_test)
  each_result =  matrix(0, nrow = n_simultaneous, ncol = 2)
  
  n_sample = n - n_test
  
  for (test in 1:n_simultaneous) {
    target_index = test
    #### select the target sample ####
    Y_test_one = Y_test[target_index]
    X_test_one = X_test[target_index, ]
    X_test_one <- matrix(X_test_one, nrow = 1)
    
    
    #### random select the one sample ####
    random_samples <- lapply(1:n_sample, function(i)
      random_index_sample(i, Y_sample, X_sample, Delta_sample))
    Y_one_time <- sapply(random_samples, function(x)
      x$y_sample)
    X_one_time <- do.call(rbind, lapply(random_samples, function(x)
      x$x_sample))
    Delta_one_time <- sapply(random_samples, function(x)
      x$delta_sample)
    
    
    ## split data to training and calibration
    n_training = floor(n_sample * 0.5)
    training_index = single_split(n_sample, n_training)
    X_training <- X_one_time[training_index$indices_part1, ]
    Y_training <- Y_one_time[training_index$indices_part1]
    Delta_training <- Delta_one_time[training_index$indices_part1]
    X_cali <- X_one_time[training_index$indices_part2, ]
    Y_cali <- Y_one_time[training_index$indices_part2]
    Delta_cali <- Delta_one_time[training_index$indices_part2]
    
    
    
    ### fit quantile function ####
    Y_obs_training <- Y_training[Delta_training == 1]
    X_obs_training <- X_training[Delta_training == 1, ]
    if (quantile_model == "grf") {
      quantile_forest <- quantile_forest(X_obs_training, Y_obs_training, quantiles = c(0.05, 0.95))
      
    } else if (quantile_model == "linear") {
      rq_model <- rq(
        Y_obs_training ~ .,
        data = as.data.frame(X_obs_training),
        tau = c(0.05, 0.95)
      )
    }
    
    
    
    ### Calculate regression on test points
    if (quantile_model == "grf") {
      quantile_test = predict(quantile_forest, X_test_one, quantiles = c(0.05, 0.95))$predictions
    } else if (quantile_model == "linear") {
      quantile_test = predict(rq_model, newdata = as.data.frame(X_test_one))
    }
    score_test = pmax(quantile_test[, 1] - y_grid, y_grid - quantile_test[, 2])
    
    
    ### Calculate regression on calibration
    Y_cali_sample_obs <- Y_cali[Delta_cali == 1]
    X_cali_sample_obs <- X_cali[Delta_cali == 1, ]
    
    if (quantile_model == "grf") {
      quantile_cali = predict(quantile_forest,
                              X_cali_sample_obs,
                              quantiles = c(0.05, 0.95))$predictions
    } else if (quantile_model == "linear") {
      quantile_cali = predict(rq_model, newdata = as.data.frame(X_cali_sample_obs))
    }
    score_cali = pmax(quantile_cali[, 1] - Y_cali_sample_obs,
                      Y_cali_sample_obs - quantile_cali[, 2])
    
    
    
    unwt_p_values = p_value(score_cali, rep(1, length(score_cali)), score_test, 1, n_grid)
    
    C_RF <- y_grid[unwt_p_values >= alpha]
    lo = min(C_RF)
    hi = max(C_RF)
    if (length(C_RF) == 0) {
      lo <- min(y_grid)
      hi <- max(y_grid)
    }
    
    in_interval <- (Y_test_one >= lo) & (Y_test_one <= hi)
    result_length = hi - lo
    
    each_result[test, ] = c(in_interval, result_length)
    
    
  }
  
  list(each_result)
  
}




stopCluster(cl)
pointwise_result = do.call(rbind, results)
res <- colMeans(pointwise_result)
cat(paste(res, collapse = " "), "\n")
