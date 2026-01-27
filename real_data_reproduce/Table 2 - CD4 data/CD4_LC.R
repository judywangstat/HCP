
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


### General parameters ###
alpha = 0.1
n_grid = 200
wt_model = "logistic"   # "logistic", "grf" 
quantile_model = "grf"   # "linear", "grf" 
n_test <- 1
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
  library(grf)
  library(quantreg)
  source(file.path(parent_dir,"my_functions.R"),chdir = TRUE)
})

results <- foreach(
  test_id = seq_len(n),  
  .combine = rbind, 
  .options.RNG = 123 
) %dorng% {

  print(test_id)
  ### Leave one out test data ###
  X_test <- X_list[[test_id]]
  Y_test <- Y_list[[test_id]]
  Delta_test <- Delta_list[[test_id]]
  n_simultaneous = length(Y_test)
  
  X_sample <- X_list[-test_id]
  Y_sample <- Y_list[-test_id]
  Delta_sample <- Delta_list[-test_id]
  
  lo_ygrid = min(unlist(Y_sample))
  up_ygrid = max(unlist(Y_sample))
  y_grid <- seq(from = lo_ygrid,
                to = up_ygrid,
                length.out = n_grid) 
  
  
  
  each_result =  matrix(0, nrow = n_simultaneous, ncol = 2)  
  n_sample = n - n_test
  for (test in 1: n_simultaneous) {
    target_index = test
    
    #### select the target sample ####
    Y_test_one = Y_test[target_index]
    X_test_one = X_test[target_index, ]
    X_test_one <- matrix(X_test_one, nrow = 1)  
    colnames(X_test_one) <- c("time", "age", "smoke", "drug", "partners", "cesd")
    
    #### random select the training data ####
    random_samples <- lapply(1:n_sample, function(i)
      random_index_sample(i, Y_sample, X_sample, Delta_sample))
    Y_one_time <- sapply(random_samples, function(x)
      x$y_sample)
    X_one_time <- do.call(rbind, lapply(random_samples, function(x)
      x$x_sample))
    Delta_one_time <- sapply(random_samples, function(x)
      x$delta_sample)
    
    
    
    
    #### Conformal Prediction Under Cross-sectional Data ####
    
    ## split data to training and calibration 
    n_sample = n - n_test
    n_training = floor(n_sample * 0.5)
    training_index = single_split(n_sample, n_training)
    X_training <- X_one_time[training_index$indices_part1,]
    Y_training <- Y_one_time[training_index$indices_part1]
    Delta_training <- Delta_one_time[training_index$indices_part1]
    X_cali <- X_one_time[training_index$indices_part2,]
    Y_cali <- Y_one_time[training_index$indices_part2]
    Delta_cali <- Delta_one_time[training_index$indices_part2]
    
    
    ### fit weight function ####
    Delta_training <- as.factor(Delta_training)
    if (wt_model == "grf") {
      prob_forest = probability_forest(X_training, Delta_training)
      
    } else if (wt_model == "logistic") {
      logistic_model <- glm(Delta_training ~ ., data = as.data.frame(X_training), family = binomial())
    }
    
    ### fit quantile function ####
    Y_obs_training <- Y_training[Delta_training == 1]
    X_obs_training <- X_training[Delta_training == 1, ]
    if (quantile_model == "grf") {
      quantile_forest <- quantile_forest(X_obs_training, Y_obs_training, quantiles = c(0.05, 0.95))
      
    } else if (quantile_model == "linear") {
      rq_model <- rq(Y_obs_training ~ ., data = as.data.frame(X_obs_training), tau=c(0.05, 0.95)) 
    }
    
    
    rm(
      Y_obs_training,
      X_obs_training, 
      Y_training,
      X_training,
      Delta_training
    )
    
    
    ### Calculate weights and quantile on test points
    if (wt_model == "grf") {
      prob_test <- predict(prob_forest, X_test_one)$predictions[, "1"]
    } else if (wt_model == "logistic") {
      prob_test <- predict(logistic_model, newdata = as.data.frame(X_test_one), type = "response")
    }
    wt_test = 1/prob_test
    max_wt <- if(missing == "20") 10 else if(missing == "50") 15 else stop("Invalid missing value")
    wt_test = pmin(wt_test, max_wt)
    
    if (quantile_model == "grf") {
      quantile_test = predict(quantile_forest, X_test_one, quantiles = c(0.05, 0.95))$predictions
    } else if (quantile_model == "linear") {
      quantile_test = predict(rq_model, newdata = as.data.frame(X_test_one))
    }
    score_test = pmax(quantile_test[,1]-y_grid, y_grid-quantile_test[,2])
    
    
    ### Calculate weights and densities on calibration
    Y_cali_sample_obs <- Y_cali[Delta_cali == 1]
    X_cali_sample_obs <- X_cali[Delta_cali == 1, ]
    
    ### Calculate weights and densities on calibration
    if (wt_model == "grf") {
      prob_cali <- predict(prob_forest, X_cali_sample_obs)$predictions[, "1"]
    } else if (wt_model == "logistic") {
      prob_cali <- predict(logistic_model, newdata = as.data.frame(X_cali_sample_obs), type = "response")
    }
    wt_cali = 1/prob_cali
    wt_cali = pmin(wt_cali, max_wt)
    
    if (quantile_model == "grf") {
      quantile_cali = predict(quantile_forest, X_cali_sample_obs, quantiles = c(0.05, 0.95))$predictions
    } else if (quantile_model == "linear") {
      quantile_cali = predict(rq_model, newdata = as.data.frame(X_cali_sample_obs))
    }
    score_cali = pmax(quantile_cali[,1] - Y_cali_sample_obs,  Y_cali_sample_obs - quantile_cali[,2])
    
    
    p_values = p_value(score_cali, wt_cali, score_test,
                       wt_test, n_grid)
    
    C_RF <- y_grid[p_values >= alpha]
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
res_sd = apply(pointwise_result, 2, sd)/sqrt(nrow(results))
cat(paste(res, collapse = " "), "\n")





