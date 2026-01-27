

# Generate random effects theta
generate_theta <- function(n, mu, sigma) {
  mvrnorm(n, mu = mu, Sigma = sigma)
}


# Generate covariate X
generate_X <- function(m, mu, sigma) {
  X_normal <- mvrnorm(m, mu = mu[-1], Sigma = sigma[-1, -1])
  #X_time <- 1:m
  X_time <- sample(1:m) / m
  X <- cbind(X_time, X_normal)
  colnames(X) <- NULL
  return(X)
}


# Generate errors epsilon
generate_epsilon <- function(m, X, scenario) {
  if (is.vector(X) && length(X) == d) {
    X <- matrix(X, nrow = 1, ncol = d)
  }
  epsilon <- vector(length = m)
  if (scenario == "Homo") {
    epsilon <- rnorm(m, 0, 1)
  } else if (scenario == "Heter") {
    epsilon <- rnorm(m, 0, 1 + 3 * abs(X[, 2]))
  } else if (scenario == "Asym") {
    epsilon <- rgamma(m, shape = 0.1, rate = 0.1)
  } else if (scenario == "Bimo" || scenario == "Bimo-fix") {
    epsilon <- ifelse(runif(m) < 0.5, rnorm(m, -8, 0.5), rnorm(m, 8, 0.5))
  }
  epsilon
}

# Generate Y
generate_Y <- function(X, theta, epsilon) {
  rowSums(X * theta) + epsilon
}

# Generate the missing data indicators
generate_missing_indicator <- function(X, beta) {
  X_tilde <- cbind(1, X)
  log_odds <- X_tilde %*% beta
  prob <- 1 / (1 + exp(-log_odds))
  rbinom(nrow(X), size = 1, prob = prob)
}


# Generate all data
generate_data <- function(n,
                          m,
                          d,
                          scenario,
                          theta_mu,
                          theta_sigma,
                          X_mu,
                          X_sigma,
                          beta) {
  theta <- generate_theta(n, theta_mu, theta_sigma)
  data_X <- vector("list", n)
  data_Y <- vector("list", n)
  data_Delta <- vector("list", n)
  for (i in 1:n) {
    X <- generate_X(m[i], X_mu, X_sigma)
    if (scenario == "Bimo-fix") {
      theta[i, ] <- rep(theta_mu[1], d)  # Fixed effects for Bimo-fix scenario
    }
    epsilon <- generate_epsilon(m[i], X, scenario)
    data_X[[i]] <- X
    data_Y[[i]] <- generate_Y(X, theta[i, ], epsilon)
    data_Delta[[i]] <- generate_missing_indicator(X, beta)
  }
  return(list(X = data_X, Y = data_Y, Delta = data_Delta))
}

# Generate single split indices for data
single_split <- function(total, count) {
  indices_part1 <- sample(1:total, size = count, replace = FALSE)
  indices_part2 <- setdiff(1:total, indices_part1)
  list(indices_part1 = indices_part1, indices_part2 = indices_part2)
}

# Generate multiple split indices for data
multi_split <- function(total, count, S_splitting) {
  splits <- replicate(S_splitting, {
    indices_part1 <- sample(1:total, size = count, replace = FALSE)
    indices_part2 <- setdiff(1:total, indices_part1)
    list(indices_part1 = indices_part1, indices_part2 = indices_part2)
  }, simplify = FALSE)
  
  indices_part1 <- lapply(splits, `[[`, "indices_part1")
  indices_part2 <- lapply(splits, `[[`, "indices_part2")
  list(indices_part1 = indices_part1, indices_part2 = indices_part2)
}



# Quotient Methods
dens_y_given_x <- function(y_val,
                           num_extra_points,
                           decay_factor,
                           y_given_x_hat,
                           hi_y_given_x,
                           lo_y_given_x,
                           h) {
  #' Used to estimate conditional density, both for linear quantile regression and quantile random forests
  #' Input:
  #' @param y_val A vector of values at which the conditional density is to be estimated.
  #' @param num_extra_points An integer specifying the number of extra points to be added to the tails for decay.
  #' @param decay_factor A numeric value specifying the decay factor for the tail density estimation.
  #' @param y_given_x_hat A matrix of predicted quantiles of the response variable `y` given `x`.
  #'                     The number of rows corresponds to the number of test/calibration samples,
  #'                     and the number of columns corresponds to the number of quantile levels.
  #' @param hi_y_given_x A matrix of upper quantile predictions of `y` given `x`. Must have the same dimensions as `y_given_x_hat`.
  #' @param lo_y_given_x A matrix of lower quantile predictions of `y` given `x`. Must have the same dimensions as `y_given_x_hat`.
  #' @param h A vector of bandwidths, one for each quantile level. The length of `h` should match the number of columns in `y_given_x_hat`.
  #'
  #' @return A matrix of estimated conditional densities. Each row corresponds to a test/calibration sample,
  #'         and each column corresponds to the estimated density at a value in `y_val`.
  
  if (!identical(dim(y_given_x_hat), dim(hi_y_given_x)) ||
      !identical(dim(hi_y_given_x), dim(lo_y_given_x))) {
    stop(
      'The matrices y_given_x_hat, hi_y_given_x, and lo_y_given_x do not have the same dimensions.'
    )
  }
  
  if (length(h) != ncol(y_given_x_hat)) {
    stop('The length of h does not match the number of columns in y_given_x_hat.')
  }
  
  # Check for crossing
  cross <- y_given_x_hat[, -1] - y_given_x_hat[, -ncol(y_given_x_hat)] <= 0
  n <- nrow(y_given_x_hat)
  cross = matrix(cross, nrow = n)
  y_given_x_hat_list <- vector("list", n)
  h_list <- vector("list", n)
  hi_y_given_x_list <- vector("list", n)
  lo_y_given_x_list <- vector("list", n)
  
  for (ii in 1:n) {
    if (any(cross[ii, ])) {
      ind <- which(cross[ii, ]) + 1
      y_given_x_hat_list[[ii]] <- y_given_x_hat[ii, -ind]
      hi_y_given_x_list[[ii]] <- hi_y_given_x[ii, -ind]
      lo_y_given_x_list[[ii]] <- lo_y_given_x[ii, -ind]
      h_list[[ii]] <- h[-ind]
    } else {
      y_given_x_hat_list[[ii]] <- y_given_x_hat[ii, ]
      hi_y_given_x_list[[ii]] <- hi_y_given_x[ii, ]
      lo_y_given_x_list[[ii]] <- lo_y_given_x[ii, ]
      h_list[[ii]] <- h
    }
  }
  
  # Density estimation
  f_val <- matrix(0, nrow = n, ncol = length(y_val))
  for (i in 1:n) {
    pdf_y_given_x <- (2 * h_list[[i]]) / (hi_y_given_x_list[[i]] - lo_y_given_x_list[[i]] - 1e-8)
    pdf_y_given_x[is.infinite(pdf_y_given_x)] <- NA
    pdf_y_given_x[pdf_y_given_x < 1e-10] <- 1e-10
    pdf_y_given_x[pdf_y_given_x > 1] <- 1
    
    valid_idx <- !is.na(pdf_y_given_x)
    x <- y_given_x_hat_list[[i]][valid_idx]
    f_x <- pdf_y_given_x[valid_idx]
    
    
    # Add points for tail decay
    gap <- quantile(diff(sort(x)), 0.5)
    if (is.na(gap)) {
      gap <- 2
    }
    #gap <- 5*max(diff(sort(x)))
    extra_points_left <- seq(min(x) - num_extra_points * gap, min(x) - gap, by = gap)
    extra_points_right <- seq(max(x) + gap, max(x) + num_extra_points * gap, by = gap)
    extra_density_left <- f_x[1] * decay_factor ^ seq(num_extra_points, 1)
    extra_density_right <- f_x[length(f_x)] * decay_factor ^ seq(1, num_extra_points)
    x_all <- c(extra_points_left, x, extra_points_right)
    f_x_all <- c(extra_density_left, f_x, extra_density_right)
    
    
    # Normalize
    normli <- f_x_all
    normli[normli < 1e-10] <- 1e-10
    normli[normli > 1] <- 1
    
    # Use unique points
    unique_x_all <- unique(x_all)
    unique_normli <- normli[match(unique_x_all, x_all)]
    
    # Filter out NA values
    valid_idx <- !is.na(unique_x_all)
    x_all_filtered <- unique_x_all[valid_idx]
    normli_filtered <- unique_normli[valid_idx]
    
    if (length(x_all_filtered) == 1 &&
        length(normli_filtered) == 1) {
      x_all_filtered <- c(x_all_filtered, x_all_filtered + 0.1)
      normli_filtered <- c(normli_filtered, normli_filtered)
    }
    
    # Interpolate
    f_val[i, ] <- approx(x_all_filtered,
                         normli_filtered,
                         y_val,
                         method = "linear",
                         rule = 2)$y
    
    # Manually adjust values outside the x range
    f_val[i, y_val < min(x_all)] <- normli[1]
    f_val[i, y_val > max(x_all)] <- normli[length(normli)]
  }
  
  return(f_val)
}

# Random subsampling
random_index_sample <- function(i, Y_cali, X_cali, Delta_cali) {
  index <- sample(seq_along(Y_cali[[i]]), 1)
  list(
    y_sample = Y_cali[[i]][index],
    x_sample = X_cali[[i]][index, , drop = FALSE],
    delta_sample = Delta_cali[[i]][index]
  )
}

random_index_sample1 <- function(i, Y_cali, X_cali, Delta_cali) {
  index <- sample(seq_along(Y_cali[[i]]), 1)
  list(
    y_sample = Y_cali[[i]][index],
    x_sample = X_cali[[i]][index, , drop = FALSE],
    delta_sample = Delta_cali[[i]][index],
    index = index
  )
}

random_index_XY <- function(i, Y_cali, X_cali) {
  if (length(Y_cali[[i]]) == 0) {
    warning(paste("Y_cali[[", i, "]] is empty"))
    return(NULL)
  }
  index <- sample(seq_along(Y_cali[[i]]), 1)
  list(y_sample = Y_cali[[i]][index], x_sample = X_cali[[i]][index, , drop = FALSE])
}

p_value <- function(R_XY, w_X, R_XY_test, w_X_test, n_grid) {
  # This function calculates p-values based on weighted comparisons between R_XY (training) and R_XY_test (test).
  # The inputs are:
  # - R_XY: A vector of residuals or scores from the calibration set.
  # - w_X: Weights associated with the calibration set.
  # - R_XY_test: A vector of residuals or scores from the test set.
  # - w_X_test: A single weight value associated with the test set.
  # - n_grid: The number of test points (should match the length of R_XY_test).
  # The function returns a vector of p-values representing the proportion of training residuals
  # greater than or equal to the test residuals, weighted by the input weights.
  
  if (length(R_XY) != length(w_X)) {
    stop('Dimension mismatch: Length of R_XY and w_X must be equal.')
  }
  if (n_grid != length(R_XY_test)) {
    stop('Dimension mismatch: Lengths of Y_grid and R_XY_test must be equal.')
  }
  
  total_w <- sum(w_X) + w_X_test
  indicators <- outer(R_XY, R_XY_test, ">=") * w_X
  p_y <- colSums(indicators) / total_w + w_X_test / total_w
  
  return(p_y)
}

get_max_wt <- function(n) {
  ifelse(n > 400, 22, ifelse(n <= 200, 2.5, 11.5))
}

quantile_levels <- function(taus, n) {
  h <- quantreg::bandwidth.rq(taus, n, TRUE, 0.1)
  if (!is.null(scenario)) {
    if (scenario == "Bimo") {
      nn <- length(h)
      portion_size <- floor(nn / 5)
      index_start <- 2 * portion_size + 1
      index_end <- 3 * portion_size
      h[index_start:index_end] = quantile(h, 0.2)
    } else if (scenario == "Bimo-fix") {
      nn <- length(h)
      portion_size <- floor(nn / 5)
      index_start <- 2 * portion_size + 1
      index_end <- 3 * portion_size
      h[index_start:index_end] = quantile(h, 0.01)
    }
  }
  taus_hi <- taus + h
  taus_lo <- taus - h
  
  for (i in seq_along(taus)) {
    while (taus_lo[i] <= 0 || taus_hi[i] >= 1) {
      h[i] <- h[i] / 2
      taus_hi[i] <- taus[i] + h[i]
      taus_lo[i] <- taus[i] - h[i]
    }
  }
  return(list(
    h = h,
    taus_hi = taus_hi,
    taus_lo = taus_lo
  ))
}



CCT_p_value_matrix <- function(p_matrix, omega) {
  # This function calculates the weighted Cauchy Combination Test (CCT) p-value for each column of a matrix.
  # The CCT is used to enhance p-values in semi-supervised settings, providing a robust method
  # to combine p-values under arbitrary dependency structures.
  #
  # Arguments:
  # - p_matrix: A numeric matrix where each column contains p-values to be combined.
  # - omega: A numeric vector of weights, with a length equal to the number of rows in p_matrix.
  #
  # The function applies the CCT formula to each column of the matrix. For each column,
  # it transforms the p-values using the Cauchy transformation, computes the weighted sum,
  # and returns a vector of combined p-values for all columns.
  #
  # Reference:
  # Liu, Yaowu, and Jun Xie. "Cauchy combination test: a powerful test with analytic p-value
  # calculation under arbitrary dependency structures." Journal of the American Statistical
  # Association 115.529 (2020): 393-402.
  #
  # Example usage:
  # p_matrix <- matrix(c(0.01, 0.03, 0.02, 0.05, 0.1, 0.07), nrow = 3, ncol = 2)
  # omega <- c(0.5, 0.3, 0.2)
  # combined_p <- CCT_p_value_matrix(p_matrix, omega)
  # print(combined_p)
  
  omega = omega[complete.cases(p_matrix)]
  
  while (any(is.na(p_matrix))) {
    p_matrix <- na.omit(p_matrix)
  }
  
  if (length(omega) != nrow(p_matrix)) {
    stop("The length of omega must match the number of rows in p_matrix.")
  }
  
  p_combined_vector <- apply(p_matrix, 2, function(p_col) {
    t0 <- sum(omega * tan((0.5 - p_col) * pi))
    p_combined <- 0.5 - (atan(t0) / pi)
    return(p_combined)
  })
  
  return(p_combined_vector)
}


Mean_p_value_matrix <- function(p_matrix) {
  # This function calculates the mean combination p-value for each column of a matrix.
  # The mean combination method provides a simple and intuitive way to combine p-values.
  #
  # Arguments:
  # - p_matrix: A numeric matrix where each column contains p-values to be combined.
  #
  # The function calculates the mean p-value for each column and returns a vector of combined p-values.
  #
  # Example usage:
  # p_matrix <- matrix(c(0.01, 0.03, 0.02, 0.05, 0.1, 0.07), nrow = 3, ncol = 2)
  # combined_p <- Mean_p_value_matrix(p_matrix)
  # print(combined_p)
  
  # Remove rows with missing values
  p_matrix <- p_matrix[complete.cases(p_matrix), ]
  
  # Calculate the mean p-value for each column
  p_combined_vector <- apply(p_matrix, 2, mean)
  
  return(p_combined_vector)
}


bound_value <- function(x, min_value, max_value) {
  pmin(pmax(x, min_value), max_value)
}


### Functions for Local P-value  ####
profile_density_lo <- function(t_grid, y_grid, cde_estimate) {
  # Check if cde_estimate is a vector or a matrix
  if (is.vector(cde_estimate)) {
    cde_estimate <- matrix(cde_estimate, nrow = 1)  # Convert to a 1-row matrix
  }
  
  # Initialize the result matrix
  g_matrix <- matrix(0, nrow = nrow(cde_estimate), ncol = length(t_grid))
  
  # Loop through each row of cde_estimate
  for (i in 1:nrow(cde_estimate)) {
    fy <- cde_estimate[i, -length(cde_estimate[i, ])]
    fy_order <- fy[order(fy)]
    Num_density <- findInterval(t_grid, fy_order)   # Num_density <- sapply(t_grid, function(t) sum(fy_order <= t))
    
    Grid_area <- (y_grid[-1] - y_grid[-length(y_grid)]) * fy
    Grid_area_order <- Grid_area[order(fy)]
    
    g <- numeric(length(t_grid))
    g[Num_density > 0] <- cumsum(Grid_area_order)[Num_density[Num_density > 0]]
    g[Num_density == 0] <- 0
    
    # Assign the result to the corresponding row in the result matrix
    g_matrix[i, ] <- g
  }
  
  return(g_matrix)
}

profile_density_up <- function(t_grid, y_grid, cde_estimate) {
  # Check if cde_estimate is a vector; if so, convert it to a matrix with 1 row
  if (is.vector(cde_estimate)) {
    cde_estimate <- matrix(cde_estimate, nrow = 1)
  }
  result <- matrix(0, nrow = nrow(cde_estimate), ncol = length(t_grid))
  
  for (i in 1:nrow(cde_estimate)) {
    fy <- cde_estimate[i, -length(cde_estimate[i, ])]
    fy_order <- fy[order(fy)]
    Num_density <- findInterval(t_grid, fy_order)   # Num_density <- sapply(t_grid, function(t) sum(fy_order <= t))
    Num_density <- length(fy_order) - Num_density
    
    Grid_area <- (y_grid[-1] - y_grid[-length(y_grid)]) * fy
    Grid_area_order <- rev(Grid_area[order(fy)])
    
    g <- numeric(length(t_grid))
    g[Num_density > 0] <- cumsum(Grid_area_order)[Num_density[Num_density > 0]]
    g[Num_density == 0] <- 0
    result[i, ] <- g
  }
  return(result)
}

which_neighbors <-  function(xTrain, xTest, k) {
  return(FNN::get.knnx(data = xTrain, query = xTest, k = k)$nn.index)
}

get_S <- function(n, scenario) {
  ifelse(n >= 300 & scenario %in% c("Bimo", "Bimo-fix"), 1, 5)
}

get_b <- function(n, scenario) {
  ifelse(n >= 300 & scenario %in% c("Bimo", "Bimo-fix"), 6, 4)
}

### Functions for Multiple Data Splitting ####
splits_K_group <- function(n, K) {
  if (K >= n) {
    stop("The number of groups K cannot be greater than the total number of data points n.")
  }
  indices <- 1:n
  shuffled_indices <- sample(indices)
  groups <- split(shuffled_indices, cut(1:n, K, labels = FALSE))
  return(groups)
}


### Calculate Conditional and Local coverage
#generate multiple Y conditional on one X, for Conditional coverage
generate_Y_given_X <- function(X, n_Y, theta_mu, theta_sigma, scenario) {
  theta_given_X <- generate_theta(n_Y, theta_mu, theta_sigma)
  epsilon_given_X <- generate_epsilon(n_Y, X, scenario)
  Y_given_X <- rowSums(sweep(theta_given_X, 2, X, `*`)) + epsilon_given_X
  return(Y_given_X)
}

find_indices <- function(X_test, x2_min, x2_max, x3_min, x3_max) {
  indices <- which(X_test$x2 >= x2_min & X_test$x2 <= x2_max &
                     X_test$x3 >= x3_min & X_test$x3 <= x3_max)
  return(indices)
}

#generate Y conditional on each X, for Local coverage: one-to-one generate
generate_Y_given_Local_X <- function(X, mu, sigma, scenario) {
  n <- nrow(X)
  theta <- generate_theta(n, mu, sigma)
  epsilon <- generate_epsilon(n, X, scenario)
  Y <- generate_Y(X, theta, epsilon)
  return(Y)
}

generate_group <- function(Ti_sample) {
  group <- rep(seq_along(Ti_sample), times = Ti_sample)
  return(group)
}