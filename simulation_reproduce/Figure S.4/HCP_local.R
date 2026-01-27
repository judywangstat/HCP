
# Load necessary libraries
library(MASS)
library(stats)
library(grf)
library(doParallel)
library(doRNG)
library(rstudioapi)
script_path <- getSourceEditorContext()$path
script_dir <- dirname(script_path)
parent_dir <- dirname(script_dir)
source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)


# General parameters
num_simulations <- 1000
n <- 300  # Number of subjects/groups
m <- rep(5, n)  # Number of observations per subject/group
d <- 4  # Dimension of X and theta
scenario = "Homo"  # "Homo", "Heter", "Asym", "Bimo", "Bimo-fix"
b_true = 4  # Number of quantile levels
alpha = 0.1
B_sampling = 1
n_grid = 200
wt_model = "logistic"   # "logistic", "grf" 
quantile_model = "linear"   # "linear", "grf" 
missing <- "20"  #c("20", "50")
#nclusters = 1


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
  X_sample <- data$X
  Y_sample <- data$Y
  Delta_sample <- data$Delta
  
  lo_ygrid = min(unlist(Y_sample))
  up_ygrid = max(unlist(Y_sample))
  y_grid <- seq(from = lo_ygrid,
                to = up_ygrid,
                length.out = n_grid)
  
  ### Given test data
  x1_values <- 0.6
  x2_values <- seq(-3, 3, length.out = 42)
  x3_values <- seq(-3, 3, length.out = 42)
  x4_values <- 0
  X_test <- expand.grid(x1 = x1_values, x2 = x2_values, x3 = x3_values, x4 = x4_values); 
  Y_test = generate_Y_given_Local_X(X_test, theta_mu, theta_sigma, scenario)
  
  
  ## Given partition
  x2_breaks <- seq(-3, 3, length.out = 7)  
  x3_breaks <- seq(-3, 3, length.out = 7)
  
  ranges <- list()
  count <- 1
  for (j in 1:(length(x3_breaks)-1)) {
    for (i in 1:(length(x2_breaks)-1)) {
      ranges[[count]] <- list(
        x2_min = x2_breaks[i], x2_max = x2_breaks[i+1],
        x3_min = x3_breaks[j], x3_max = x3_breaks[j+1]
      )
      count <- count + 1
    }
  }
  
  indices_list <- lapply(ranges, function(r) {
    find_indices(X_test, r$x2_min, r$x2_max, r$x3_min, r$x3_max)
  })
  colnames(X_test) <- NULL; X_test = as.matrix(X_test)
  
  
  
  ###   split data sample to training and calibration ####
  n_training = floor(n * 0.3)
  training_index = single_split(n, n_training)
  X_training <- X_sample[training_index$indices_part1]
  Y_training <- Y_sample[training_index$indices_part1]
  Delta_training <- Delta_sample[training_index$indices_part1]
  X_cali <- X_sample[training_index$indices_part2]
  Y_cali <- Y_sample[training_index$indices_part2]
  Delta_cali <- Delta_sample[training_index$indices_part2]
  
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
    logistic_model <- glm(Delta_whole_training ~ ., data = as.data.frame(X_whole_training), family = binomial())
  }
  
  
  
  ### fit conditional density function ####
  taus <- (1:(2 ^ b_true - 1)) / (2 ^ b_true)
  taus_h = quantile_levels(taus, length(Y_obs_training))
  h = taus_h$h
  taus_hi <- taus_h$taus_hi
  taus_lo <- taus_h$taus_lo
  all_taus <- sort(c(taus_lo, taus, taus_hi))
  
  if (quantile_model == "grf") {
    quantile_forest <- quantile_forest(X_obs_training, Y_obs_training)
  } else if (quantile_model == "linear") {
    rq_model <- rq(Y_obs_training ~ ., data = as.data.frame(X_obs_training), tau=all_taus) 
  }
  
  
  
  #####  fit space partition  #####
  ### calculate density on training data, used for kmeans
  if (nclusters > 1) {
    if (quantile_model == "grf") {
      YFit_LQR_train = predict(quantile_forest, X_whole_training, quantiles = taus)$predictions
      YFit_hi_LQR_train = predict(quantile_forest, X_whole_training, quantiles = taus_hi)$predictions
      YFit_lo_LQR_train = predict(quantile_forest, X_whole_training, quantiles = taus_lo)$predictions
    } else if (quantile_model == "linear") {
      predict_rq_train = predict(rq_model, newdata = as.data.frame(X_whole_training))
      YFit_LQR_train = predict_rq_train[, match(taus, all_taus)]
      YFit_hi_LQR_train = predict_rq_train[, match(taus_hi, all_taus)]
      YFit_lo_LQR_train = predict_rq_train[, match(taus_lo, all_taus)]
    }
    dens_Pyx_LQR_train = dens_y_given_x(y_grid,
                                        20,
                                        0.7,
                                        YFit_LQR_train,
                                        YFit_hi_LQR_train,
                                        YFit_lo_LQR_train,
                                        h)
    
    ### calculate k-means center based on profile distance 
    t_num = 500
    t_grid <- seq(0, max(dens_Pyx_LQR_train), length.out = t_num)
    profile_train = profile_density_lo(t_grid, y_grid, dens_Pyx_LQR_train)
    kmeans_result <- kmeans(profile_train, centers = nclusters)
    kmeans_centers <- kmeans_result$centers
    
    rm(
      data,
      training_index,
      X_sample,
      Y_sample,
      Delta_sample,
      X_whole_training,
      Y_whole_training,
      Y_obs_training,
      X_obs_training
    )
    
  }
  
  ### Calculate weights and density on test points
  if (wt_model == "grf") {
    prob_test <- predict(prob_forest, X_test)$predictions[, "1"]
  } else if (wt_model == "logistic") {
    prob_test <- predict(logistic_model, newdata = as.data.frame(X_test), type = "response")
  }
  wt_test = 1/prob_test
  wt_test = pmin(wt_test, max_wt)
  
  
  
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
  
  ### Find X_{n+1} in which partition 
  if (nclusters > 1) {
    profile_test = profile_density_lo(t_grid, y_grid, dens_Pyx_LQR_te)
    which_partition_test = which_neighbors(kmeans_centers, profile_test, 1)
    
    rm(YFit_LQR_te, YFit_hi_LQR_te, YFit_lo_LQR_te)
  }
  
  ### Repeat subsampling the calibration data  ####
  n_cali = n - n_training
  n_simultaneous = length(Y_test)
  p_values <- vector("list", n_simultaneous)
  for (ii in 1:n_simultaneous) {
    p_values[[ii]] <- matrix(NA, nrow = B_sampling, ncol = n_grid)  
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
    
    rm(random_samples,
       Y_cali_sample,
       X_cali_sample,
       Delta_cali_sample)
    
    
    
    #### Find calibration data which in partition A(x_{n+1})  #### 
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
    if (nclusters > 1) {
    dens_Pyx_LQR_ca_grid = dens_y_given_x(y_grid,
                                          20,
                                          0.7,
                                          YFit_LQR_ca,
                                          YFit_hi_LQR_ca,
                                          YFit_lo_LQR_ca,
                                          h)
    profile_cali = profile_density_lo(t_grid, y_grid, dens_Pyx_LQR_ca_grid)
    which_partition_cali = which_neighbors(kmeans_centers, profile_cali, 1)
    
    index_A_x = sapply(which_partition_test, function(x) which_partition_cali %in% x)
    
    rm(
      YFit_LQR_ca,
      YFit_hi_LQR_ca,
      YFit_lo_LQR_ca,
      dens_Pyx_LQR_ca_grid
    )
    } else if (nclusters == 1) {
      dens_Pyx_LQR_ca_full = dens_y_given_x(Y_cali_sample_obs,
                                            20,
                                            0.7,
                                            YFit_LQR_ca,
                                            YFit_hi_LQR_ca,
                                            YFit_lo_LQR_ca,
                                            h)
      dens_Pyx_LQR_ca = diag(dens_Pyx_LQR_ca_full)
      
      ### Calculate weights
      if (wt_model == "grf") {
        prob_cali <- predict(prob_forest, X_cali_sample_obs)$predictions[, "1"]
      } else if (wt_model == "logistic") {
        prob_cali <- predict(logistic_model,
                             newdata = as.data.frame(X_cali_sample_obs),
                             type = "response")
      }
      wt_cali = 1 / prob_cali
      wt_cali = pmin(wt_cali, max_wt)
    }
    
    
    
    
    
    
    for (same in 1: n_simultaneous) {
      if (nclusters > 1) {
      ## select local sample
      if (length(which(index_A_x[,same] == 1)) <= 1) {
        warning("Warning: partition A(x_{n+1}) is empty. Skipping iteration.")
        next 
      }  
      Y_cali_sample_obs_loc = Y_cali_sample_obs[index_A_x[,same]]
      X_cali_sample_obs_loc = X_cali_sample_obs[index_A_x[,same], ] 
      
      
      ### Calculate weights and densities on Local_calibration
      if (wt_model == "grf") {
        prob_cali <- predict(prob_forest, X_cali_sample_obs_loc)$predictions[, "1"]
      } else if (wt_model == "logistic") {
        prob_cali <- predict(logistic_model, newdata = as.data.frame(X_cali_sample_obs_loc), type = "response")
      }
      wt_cali = 1/prob_cali
      wt_cali = pmin(wt_cali, max_wt)
      
      if (quantile_model == "grf") {
        YFit_LQR_ca_loc = predict(quantile_forest, X_cali_sample_obs_loc, quantiles = taus)$predictions
        YFit_hi_LQR_ca_loc = predict(quantile_forest, X_cali_sample_obs_loc, quantiles = taus_hi)$predictions
        YFit_lo_LQR_ca_loc = predict(quantile_forest, X_cali_sample_obs_loc, quantiles = taus_lo)$predictions
      } else if (quantile_model == "linear") {
        predict_rq_ca = predict(rq_model, newdata = as.data.frame(X_cali_sample_obs_loc))
        YFit_LQR_ca_loc = predict_rq_ca[, match(taus, all_taus)]
        YFit_hi_LQR_ca_loc = predict_rq_ca[, match(taus_hi, all_taus)]
        YFit_lo_LQR_ca_loc = predict_rq_ca[, match(taus_lo, all_taus)]
      }
      dens_Pyx_LQR_ca_full = dens_y_given_x(Y_cali_sample_obs_loc,
                                            20,
                                            0.7,
                                            YFit_LQR_ca_loc,
                                            YFit_hi_LQR_ca_loc,
                                            YFit_lo_LQR_ca_loc,
                                            h)
      dens_Pyx_LQR_ca = diag(dens_Pyx_LQR_ca_full)
      
      
      
      ## p-value at each new point
      p_values[[same]][b_sampling, ] = p_value(-dens_Pyx_LQR_ca,
                                               wt_cali,-dens_Pyx_LQR_te[same,],
                                               wt_test[same], n_grid)
      } else if (nclusters == 1) {
        p_values[[same]][b_sampling, ] = p_value(-dens_Pyx_LQR_ca,
                                                 wt_cali,
                                                 -dens_Pyx_LQR_te[same, ],
                                                 wt_test[same],
                                                 n_grid)
        
      }
      
    }
  }
  
  each_result =  matrix(0, nrow = n_simultaneous, ncol = 2)  
  omega = rep(1 / B_sampling, B_sampling)
  for (same in 1: n_simultaneous) {
    p_combined <- CCT_p_value_matrix(p_values[[same]], omega)
    C_RF <- y_grid[p_combined >= alpha]
    lo = min(C_RF)
    hi = max(C_RF)
    if (length(C_RF) == 0) {
      lo <- min(y_grid)
      hi <- max(y_grid)
    }
    in_interval <- (Y_test[same] >= lo) & (Y_test[same] <= hi)
    result_length = hi - lo
    if (scenario == "Bimo" || scenario == "Bimo-fix" ) {
      if(length(which(diff(C_RF) > 1))>0){
        differences = which(diff(C_RF) > 1)[1]
        lo_dif = C_RF[differences]
        hi_dif =C_RF[differences+1]
        leng_dif = hi_dif - lo_dif
        result_length = result_length - leng_dif
        
        in_interval <- (Y_test[same] >= lo & Y_test[same] <= lo_dif) || (Y_test[same] >= hi_dif & Y_test[same] <= hi)
      }
    }
    
    each_result[same, ] = c(in_interval, result_length)
  }
  
  ### calculate local coverage
  local_results <- lapply(indices_list, function(indices) {
    each_result[indices, , drop = FALSE]  
  })
  
  local_means <- lapply(local_results, function(matrix) {
    colMeans(matrix, na.rm = TRUE)  
  })
  
  local_means = unlist(local_means)
  
  #results[[sim]] = each_result
  list(local_means)
  
  
}


pointwise_result = do.call(rbind, results) 
mean_result = colMeans(pointwise_result)
coverage_values <- mean_result[seq(1, ncol(pointwise_result), by = 2)]
length_values <- mean_result[-seq(1, ncol(pointwise_result), by = 2)]

final_matrix <- rbind(coverage_values, length_values)
rownames(final_matrix) <- c("Coverage", "Length")
print(final_matrix)






