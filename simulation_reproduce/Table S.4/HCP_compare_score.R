

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
#n <- 300  # Number of subjects/groups
n_test <- 1
m <- rep(5, n)  # Number of observations per subject/group
d <- 4  # Dimension of X and theta
#scenario = "Asym"  # "Homo", "Heter", "Asym", "Bimo", "Bimo-fix"
b_true = 4  # Number of quantile levels
B_sampling = 5
S_splitting = 5
alpha = 0.1
n_grid = 200
wt_model = "logistic"   # "logistic", "grf"
quantile_model = "linear"   # "linear", "grf"
missing <- "20"  #c("20", "50")
reg_model = "linear"   # "linear", "grf" 


# Parameters for distributions
theta_mu <- rep(2, d)
theta_sigma <- diag(1, d)
X_mu <- rep(0, d)
X_sigma <- diag(1, d)
beta <- c(3, 0, 2, 2, 2)
max_wt <- 30


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
  
  
  
  ###  multi_split data sample to training and calibration ####
  n_sample = n - n_test
  n_training = floor(n_sample * 0.5)
  S_splitting = get_S(n, scenario)
  multi_index = multi_split(n_sample, n_training, S_splitting)
  train_multi_index = multi_index$indices_part1
  cali_multi_index = multi_index$indices_part2
  
  if (!(
    length(train_multi_index) == length(cali_multi_index) &&
    length(cali_multi_index) == S_splitting
  )) {
    stop(
      "length(train_multi_index), length(cali_multi_index), and S_splitting are not equal. Program has been stopped."
    )
  }
  
  p_values_S_group <- vector("list", S_splitting)
  res_p_values_K_group <- vector("list", S_splitting)
  ## create p-value matrix:  1. test point 2. B_sampling 3. y_grid
  for (ii in 1:S_splitting) {
    p_values_S_group[[ii]] <- matrix(0, nrow = n_simultaneous, ncol = n_grid)
    res_p_values_K_group[[ii]] <- matrix(0, nrow = n_simultaneous, ncol = n_grid)
  }
  
  
  for (s_split in 1:S_splitting) {
    ### select data to two part: training and calibration ####
    training_index = train_multi_index[[s_split]]
    calibration_index = cali_multi_index[[s_split]]
    
    X_training <- X_sample[training_index]
    Y_training <- Y_sample[training_index]
    Delta_training <- Delta_sample[training_index]
    n_training = length(training_index)
    
    
    X_cali <- X_sample[calibration_index]
    Y_cali <- Y_sample[calibration_index]
    Delta_cali <- Delta_sample[calibration_index]
    n_cali = length(calibration_index)
    
    Y_whole_training <- unlist(Y_training)
    Delta_whole_training <- unlist(Delta_training)
    X_whole_training <- do.call(rbind, X_training)
    Y_obs_training <- Y_whole_training[Delta_whole_training == 1]
    X_obs_training <- X_whole_training[Delta_whole_training == 1, ]
    
    
    ### fit weight function ####
    Delta_whole_training <- as.factor(Delta_whole_training)
    if (wt_model == "grf") {
      prob_forest = probability_forest(X_whole_training, Delta_whole_training)
      
    } else if (wt_model == "logistic") {
      logistic_model <- glm(
        Delta_whole_training ~ .,
        data = as.data.frame(X_whole_training),
        family = binomial()
      )
    }
    
    
    ### fit conditional density function ####
    b_true = get_b(n, scenario) 
    taus <- (1:(2^b_true - 1)) / (2^b_true)
    taus_h = quantile_levels(taus, length(Y_obs_training))
    h = taus_h$h
    taus_hi <- taus_h$taus_hi
    taus_lo <- taus_h$taus_lo
    all_taus <- sort(c(taus_lo, taus, taus_hi))
    
    if (quantile_model == "grf") {
      quantile_forest <- quantile_forest(X_obs_training, Y_obs_training)
      #quantile_forest <- quantile_forest(X_obs_training, Y_obs_training, quantiles = all_taus)
      
    } else if (quantile_model == "linear") {
      rq_model <- rq(Y_obs_training ~ .,
                     data = as.data.frame(X_obs_training),
                     tau = all_taus)
    }
    
    
    ### fit regression function ####
    if (reg_model == "grf") {
      regression_forest <- regression_forest(X_obs_training, Y_obs_training)
    } else if (reg_model == "linear") {
      lm_model <- lm(Y_obs_training ~ ., data = as.data.frame(X_obs_training))
    }
    
    
    ### Calculate weights and density on test points
    if (wt_model == "grf") {
      prob_test <- predict(prob_forest, X_test)$predictions[, "1"]
    } else if (wt_model == "logistic") {
      prob_test <- predict(logistic_model,
                           newdata = as.data.frame(X_test),
                           type = "response")
    }
    wt_test = 1 / prob_test
    wt_test = pmin(wt_test, max_wt)
    
    ### Calculate densities on test points
    if (quantile_model == "grf") {
      YFit_LQR_te = predict(quantile_forest, X_test, quantiles = taus)$predictions
      YFit_hi_LQR_te = predict(quantile_forest, X_test, quantiles = taus_hi)$predictions
      YFit_lo_LQR_te = predict(quantile_forest, X_test, quantiles = taus_lo)$predictions
    } else if (quantile_model == "linear") {
      predict_rq_test = predict(rq_model, newdata = as.data.frame(X_test))
      YFit_LQR_te = predict_rq_test[, match(taus, all_taus)]
      YFit_hi_LQR_te = predict_rq_test[, match(taus_hi, all_taus)]
      YFit_lo_LQR_te = predict_rq_test[, match(taus_lo, all_taus)]
    }
    dens_Pyx_LQR_te = dens_y_given_x(y_grid,
                                     20,
                                     0.7,
                                     YFit_LQR_te,
                                     YFit_hi_LQR_te,
                                     YFit_lo_LQR_te,
                                     h)
    
    ### Calculate residuals on test points
    if (reg_model == "grf") {
      reg_test = predict(regression_forest, X_test)$predictions
    } else if (reg_model == "linear") {
      reg_test <- predict(lm_model, newdata = as.data.frame(X_test))
    }
    score_test = abs(outer(reg_test, y_grid, "-"))
    
    
    
    
    
    ### Repeat subsampling the calibration data
    p_values <- vector("list", n_simultaneous)
    res_p_values <- vector("list", n_simultaneous)
    
    ## create p-value matrix:  1. test point 2. B_sampling 3. y_grid
    for (ii in 1:n_simultaneous) {
      p_values[[ii]] <- matrix(0, nrow = B_sampling, ncol = n_grid)
      res_p_values[[ii]] <- matrix(0, nrow = B_sampling, ncol = n_grid)
    }
    
    
    for (b_sampling in 1:B_sampling) {
      #first subsampling, then using observed data
      random_samples <- lapply(1:n_cali, function(i)
        random_index_sample(i, Y_cali, X_cali, Delta_cali))
      Y_cali_sample <- sapply(random_samples, function(x)
        x$y_sample)
      X_cali_sample <- do.call(rbind, lapply(random_samples, function(x)
        x$x_sample))
      Delta_cali_sample <- sapply(random_samples, function(x)
        x$delta_sample)
      
      Y_cali_sample_obs <- Y_cali_sample[Delta_cali_sample == 1]
      X_cali_sample_obs <- X_cali_sample[Delta_cali_sample == 1, ]
      n_cali_obs = length(Y_cali_sample_obs)
      
      rm(random_samples,
         Y_cali_sample,
         X_cali_sample,
         Delta_cali_sample)
      
      ### Calculate weights and densities on calibration
      if (wt_model == "grf") {
        prob_cali <- predict(prob_forest, X_cali_sample_obs)$predictions[, "1"]
      } else if (wt_model == "logistic") {
        prob_cali <- predict(logistic_model,
                             newdata = as.data.frame(X_cali_sample_obs),
                             type = "response")
      }
      wt_cali = 1 / prob_cali
      wt_cali = pmin(wt_cali, max_wt)
      
      if (quantile_model == "grf") {
        YFit_LQR_ca = predict(quantile_forest, X_cali_sample_obs, quantiles = taus)$predictions
        YFit_hi_LQR_ca = predict(quantile_forest, X_cali_sample_obs, quantiles = taus_hi)$predictions
        YFit_lo_LQR_ca = predict(quantile_forest, X_cali_sample_obs, quantiles = taus_lo)$predictions
      } else if (quantile_model == "linear") {
        predict_rq_ca = predict(rq_model, newdata = as.data.frame(X_cali_sample_obs))
        YFit_LQR_ca = predict_rq_ca[, match(taus, all_taus)]
        YFit_hi_LQR_ca = predict_rq_ca[, match(taus_hi, all_taus)]
        YFit_lo_LQR_ca = predict_rq_ca[, match(taus_lo, all_taus)]
      }
      dens_Pyx_LQR_ca_full = dens_y_given_x(Y_cali_sample_obs,
                                            20,
                                            0.7,
                                            YFit_LQR_ca,
                                            YFit_hi_LQR_ca,
                                            YFit_lo_LQR_ca,
                                            h)
      dens_Pyx_LQR_ca = diag(dens_Pyx_LQR_ca_full)
      
      
      ### Calculate residuals on calibration
      if (reg_model == "grf") {
        reg_cali = predict(regression_forest, X_cali_sample_obs)$predictions
      } else if (reg_model == "linear") {
        reg_cali <- predict(lm_model, newdata = as.data.frame(X_cali_sample_obs))
      }
      score_cali = abs(Y_cali_sample_obs - reg_cali)
      
      
      
      ## p-value at each new point
      for (same in 1:n_simultaneous) {
        p_values[[same]][b_sampling, ] = p_value(-dens_Pyx_LQR_ca,
                                                 wt_cali,-dens_Pyx_LQR_te[same, ],
                                                 wt_test[same],
                                                 n_grid)
        res_p_values[[same]][b_sampling, ] = p_value(score_cali,
                                                     wt_cali, score_test[same,],
                                                     wt_test[same], n_grid)
      }
    }
    
    
    ## Combine B-sampling p-values
    omega_B_sampling = rep(1 / B_sampling, B_sampling)
    for (same in 1:n_simultaneous) {
      p_values_S_group[[s_split]][same, ] = CCT_p_value_matrix(p_values[[same]], omega_B_sampling)
      res_p_values_K_group[[s_split]][same,] = CCT_p_value_matrix(res_p_values[[same]], omega_B_sampling)
    }
  }
  
  
  each_result =  matrix(0, nrow = n_simultaneous, ncol = 4)
  omega_S_group = rep(1 / S_splitting, S_splitting)
  for (same in 1:n_simultaneous) {
    ## Combine K-group p-values
    p_S_group = do.call(rbind, lapply(p_values_S_group, function(x)
      x[same, ]))
    p_combined <- CCT_p_value_matrix(p_S_group, omega_S_group)
    C_RF <- y_grid[p_combined >= alpha]
    lo = min(C_RF)
    hi = max(C_RF)
    if (length(C_RF) == 0) {
      lo <- min(y_grid)
      hi <- max(y_grid)
    }
    in_interval <- (Y_test[same] >= lo) & (Y_test[same] <= hi)
    result_length = hi - lo
    if (scenario == "Bimo" || scenario == "Bimo-fix") {
      if (length(which(diff(C_RF) > 1)) > 0) {
        differences = which(diff(C_RF) > 1)[1]
        lo_dif = C_RF[differences]
        hi_dif = C_RF[differences + 1]
        leng_dif = hi_dif - lo_dif
        result_length = result_length - leng_dif
        
        in_interval <- (Y_test[same] >= lo &
                          Y_test[same] <= lo_dif) ||
          (Y_test[same] >= hi_dif & Y_test[same] <= hi)
      }
    }
    
    res_p_K_group = do.call(rbind, lapply(res_p_values_K_group, function(x) x[same, ]))  
    res_p_combined <- CCT_p_value_matrix( res_p_K_group, omega_S_group)
    res_C_RF <- y_grid[res_p_combined >= alpha]
    res_lo = min(res_C_RF)
    res_hi = max(res_C_RF)
    if (length(C_RF) == 0) {
      res_lo <- min(y_grid)
      res_hi <- max(y_grid)
    }
    res_in_interval <- (Y_test[same] >= res_lo) & (Y_test[same] <= res_hi)
    res_result_length = res_hi - res_lo
    
    
    each_result[same, ] = c(in_interval, result_length, res_in_interval, res_result_length)
  }
  
  
  #results[[sim]] = each_result
  list(each_result)
}


stopCluster(cl)

pointwise_result = do.call(rbind, results)
res <- colMeans(pointwise_result)
cat(paste(res, collapse = " "), "\n")

