# Load necessary libraries
library(MASS)
library(stats)
library(grf)
library(doParallel)
library(doRNG)
library(quantreg)
library(grf)
library(randomForest)

library(rstudioapi)
script_path <- getSourceEditorContext()$path
script_dir <- dirname(script_path)
parent_dir <- dirname(script_dir)
source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)


### General parameters ###
alpha = 0.1
n_grid = 200
n_test <- 1
B_sampling = 1
reg_model = "linear"     # "linear", "grf" 

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
  ### Repeat subsampling to a cross sectional data####
  n_sample = n - n_test
  unwt_p_values <- vector("list", n_simultaneous)
  for (ii in 1:n_simultaneous) {
    unwt_p_values[[ii]] <- matrix(NA, nrow = B_sampling, ncol = n_grid)  
  }
  
  for (b_sampling in 1:B_sampling) {
    random_samples <- lapply(1:n_sample, function(i)
      random_index_sample(i, Y_sample, X_sample, Delta_sample))
    Y_sample_one <- sapply(random_samples, function(x)
      x$y_sample)
    X_sample_one <- do.call(rbind, lapply(random_samples, function(x)
      x$x_sample))
    Delta_sample_one <- sapply(random_samples, function(x)
      x$delta_sample)
    
    
    ### split data to training and calibration ####
    n_training = floor(n_sample * 0.5)
    training_index = single_split(n_sample, n_training)
    X_training <- X_sample_one[training_index$indices_part1,]
    Y_training <- Y_sample_one[training_index$indices_part1]
    Delta_training <- Delta_sample_one[training_index$indices_part1]
    X_cali <- X_sample_one[training_index$indices_part2,]
    Y_cali <- Y_sample_one[training_index$indices_part2]
    Delta_cali <- Delta_sample_one[training_index$indices_part2]
    
    ### fit regression function ####
    Y_obs_training <- Y_training[Delta_training == 1]
    X_obs_training <- X_training[Delta_training == 1, ]
    
    if (reg_model == "grf") {
      regression_forest <- regression_forest(X_obs_training, Y_obs_training)
    } else if (reg_model == "linear") {
      lm_model <- lm(Y_obs_training ~ ., data = as.data.frame(X_obs_training))
    }
    
    rm(
      random_samples,
      training_index,
      X_obs_training, 
      Y_training,
      X_training,
      Delta_training
    )
    
    ### Calculate regression on test points
    if (reg_model == "grf") {
      reg_test = predict(regression_forest, X_test)$predictions
    } else if (reg_model == "linear") {
      reg_test <- predict(lm_model, newdata = as.data.frame(X_test))
    }
    score_test = abs(outer(reg_test, y_grid, "-"))
    
    ### Calculate regression on calibration
    Y_cali_sample_obs <- Y_cali[Delta_cali == 1]
    X_cali_sample_obs <- X_cali[Delta_cali == 1, ]
    n_cali_obs = length(Y_cali_sample_obs)
    
    
    if (reg_model == "grf") {
      reg_cali = predict(regression_forest, X_cali_sample_obs)$predictions
    } else if (reg_model == "linear") {
      reg_cali <- predict(lm_model, newdata = as.data.frame(X_cali_sample_obs))
    }
    score_cali = abs(Y_cali_sample_obs - reg_cali)
    
    
    
    
    ## p-value
    for (same in 1: n_simultaneous) {
      unwt_p_values[[same]][b_sampling, ] = p_value(score_cali,
                                                    rep(1, n_cali_obs) , score_test[same,],
                                                    1, n_grid)
      # plot(y_grid, omni_p_values[[same]][b_sampling, ])
      # lines(y_grid, unwt_p_values[[same]][b_sampling, ], col = "red")
      # abline(h=0.1, col = "blue")
    }
    
    rm(Y_cali,
       X_cali,
       Delta_cali,
       Y_cali_sample_obs,
       X_cali_sample_obs
    )
  }
  
  
  each_result =  matrix(0, nrow = n_simultaneous, ncol = 2)  
  for (same in 1: n_simultaneous) {
    unwt_p_combined <- colMeans(unwt_p_values[[same]], na.rm = TRUE)
    unwt_C_RF <- y_grid[unwt_p_combined >= alpha]
    unwt_lo = min(unwt_C_RF)
    unwt_hi = max(unwt_C_RF)
    if (length(unwt_C_RF) == 0) {
      unwt_lo <- min(y_grid)
      unwt_hi <- max(y_grid)
    }
    unwt_in_interval <- (Y_test[same] >= unwt_lo) & (Y_test[same] <= unwt_hi)
    unwt_result_length = unwt_hi - unwt_lo
    
    
    each_result[same, ] = c(unwt_in_interval, unwt_result_length)
  }
  
  
  
  list(each_result)
}

stopCluster(cl)


pointwise_result = do.call(rbind, results)
res <- colMeans(pointwise_result)
res_sd = apply(pointwise_result, 2, sd)/sqrt(nrow(results))
cat(paste(res, collapse = " "), "\n")






