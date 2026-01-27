# Load necessary libraries
library(MASS)
library(stats)
library(grf)
library(doParallel)
library(doRNG)
library(quantreg)
library(randomForest)

library(rstudioapi)
script_path <- getSourceEditorContext()$path
script_dir <- dirname(script_path)
parent_dir <- dirname(script_dir)
source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)


# General parameters
scenario = NULL
n_test <- 1
b_true = 4  # Number of quantile levels
alpha = 0.1/4
B_sampling = 5
S_splitting = 5
n_grid = 200
wt_model = "logistic"   # "logistic", "grf" 
quantile_model = "linear"   # "linear", "grf" 



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
Y_pred <- m_X_pred +  1.2 * r_X_pred * e
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
  library(quantreg)
  source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)
})

results <- foreach(
  test_id = seq_len(n),  
  .combine = rbind,  
  .options.RNG = 123  
) %dorng% {
  ### Leave one out test data ###
  print(test_id)
  X_test <- X_list[[test_id]]
  Y_test <- Y_list_imputed[[test_id]]
  Delta_test <- Delta_list[[test_id]]
  n_simultaneous = length(Y_test)
  
  X_sample <- X_list[-test_id]
  Y_sample <- Y_list[-test_id]
  Delta_sample <- Delta_list[-test_id]
  
  
  lo_ygrid = min(unlist(Y_sample), na.rm = TRUE)
  up_ygrid = max(unlist(Y_sample), na.rm = TRUE)
  y_grid <- seq(from = lo_ygrid,
                to = up_ygrid,
                length.out = n_grid) 
  
  
  ###  multi_split data sample to training and calibration ####
  n_sample = n - n_test
  n_training = floor(n_sample * 0.5)
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
  ## create p-value matrix:  1. test point 2. B_sampling 3. y_grid
  for (ii in 1:S_splitting) {
    p_values_S_group[[ii]] <- matrix(0, nrow = n_simultaneous, ncol = n_grid)
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
    
    
    
    ### Calculate weights and density on test points
    if (wt_model == "grf") {
      prob_test <- predict(prob_forest, X_test)$predictions[, "1"]
    } else if (wt_model == "logistic") {
      prob_test <- predict(logistic_model,
                           newdata = as.data.frame(X_test),
                           type = "response")
    }
    wt_test = 1 / prob_test
    wt_test = pmin(wt_test, 50)
    
    
    if (quantile_model == "grf") {
      YFit_LQR_te = predict(quantile_forest, X_test, quantiles = taus)$predictions
      YFit_hi_LQR_te = predict(quantile_forest, X_test, quantiles = taus_hi)$predictions
      YFit_lo_LQR_te = predict(quantile_forest, X_test, quantiles = taus_lo)$predictions
    } else if (quantile_model == "linear") {
      predict_rq_test = predict(rq_model, newdata = as.data.frame(X_test))
      YFit_LQR_te = predict_rq_test[, match(taus, all_taus)]; YFit_LQR_te = matrix(YFit_LQR_te, nrow = n_simultaneous)
      YFit_hi_LQR_te = predict_rq_test[, match(taus_hi, all_taus)]; YFit_hi_LQR_te = matrix(YFit_hi_LQR_te, nrow = n_simultaneous)
      YFit_lo_LQR_te = predict_rq_test[, match(taus_lo, all_taus)]; YFit_lo_LQR_te = matrix(YFit_lo_LQR_te, nrow = n_simultaneous)
    }
    dens_Pyx_LQR_te = dens_y_given_x(y_grid,
                                     20,
                                     0.7,
                                     YFit_LQR_te,
                                     YFit_hi_LQR_te,
                                     YFit_lo_LQR_te,
                                     h)
    
    
    ### Repeat subsampling the calibration data
    p_values <- vector("list", n_simultaneous)
    unwt_p_values <- vector("list", n_simultaneous)
    
    ## create p-value matrix:  1. test point 2. B_sampling 3. y_grid
    for (ii in 1:n_simultaneous) {
      p_values[[ii]] <- matrix(0, nrow = B_sampling, ncol = n_grid)
      unwt_p_values[[ii]] <- matrix(0, nrow = B_sampling, ncol = n_grid)
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
      #prob_cali  = bound_value(prob_cali, 0.05, 0.95)
      wt_cali = 1 / prob_cali
      wt_cali = pmin(wt_cali, 50)
      
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
      
      
      
      ## p-value at each new point
      for (same in 1:n_simultaneous) {
        p_values[[same]][b_sampling, ] = p_value(-dens_Pyx_LQR_ca,
                                                 wt_cali,-dens_Pyx_LQR_te[same, ],
                                                 wt_test[same],
                                                 n_grid)
      }
    }
    
    
    ## Combine B-sampling p-values
    omega_B_sampling = rep(1 / B_sampling, B_sampling)
    for (same in 1:n_simultaneous) {
      p_values_S_group[[s_split]][same, ] = CCT_p_value_matrix(p_values[[same]], omega_B_sampling)
    }
  }
  
  
  each_result =  matrix(0, nrow = n_simultaneous, ncol = 2)
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
    each_result[same, ] = c(in_interval, result_length)
  }
  
  
  #results[[test_id]] = each_result
  list(each_result)
}

stopCluster(cl)
check_results = matrix(0, nrow = n, ncol = 2)  
check_results[,1] = sapply(results, function(mat) {
  if (all(mat[, 1] == 1)) {
    return(1)  
  } else {
    return(0) 
  }
})

check_results[,2] = sapply(results, function(mat) {
  return(mean(mat[, 2]))  
})



res <- colMeans(check_results)
cat(paste(res, collapse = " "), "\n")

