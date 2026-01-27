# ==============================================
# Oracle simulation: average coverage and length
# ==============================================

library(MASS)
library(stats)
library(doParallel)
library(doRNG)
library(lme4)
library(merTools)
library(rstudioapi)

script_path <- getSourceEditorContext()$path
script_dir  <- dirname(script_path)
parent_dir  <- dirname(script_dir)
source(file.path(parent_dir,"my_functions.R"), chdir = TRUE)

# ----------------------------------------------
# Oracle percentile intervals
# ----------------------------------------------
oracle_pop_interval <- function(X,
                                scenario = c("Homo","Heter","Asym","Bimo","Bimo-fix"),
                                alpha = 0.10,
                                B = 1e5,
                                seed = NULL) {
  scenario <- match.arg(scenario)
  if (is.null(dim(X))) { stopifnot(length(X) == 4); X <- matrix(X, nrow = 1) }
  stopifnot(ncol(X) == 4)
  n <- nrow(X)
  if (!is.null(seed)) set.seed(seed)
  
  mu <- as.numeric(X %*% rep(2, 4))
  v_theta <- rowSums(X * X)
  
  # Normal-based interval
  norm_int <- function(mu, var, alpha) {
    zL <- qnorm(alpha/2); zU <- qnorm(1 - alpha/2)
    sd <- sqrt(var)
    cbind(mu + sd * zL, mu + sd * zU)
  }
  
  # Homoscedastic normal errors
  if (scenario == "Homo") {
    out <- norm_int(mu, v_theta + 1, alpha)
    colnames(out) <- c("Lower_PI","Upper_PI")
    return(out)
  }
  
  # Heteroscedastic normal errors
  if (scenario == "Heter") {
    var_eps <- as.numeric((1 + 3 * abs(X[, 2]))^2)
    out <- norm_int(mu, v_theta + var_eps, alpha)
    colnames(out) <- c("Lower_PI","Upper_PI")
    return(out)
  }
  
  # Fixed bimodal mixture
  if (scenario == "Bimo-fix") {
    comp <- sample.int(2, size = B, replace = TRUE)
    means <- ifelse(comp == 1, -8, 8)
    eps <- rnorm(B, mean = means, sd = 0.5)
    Q <- t(vapply(mu, function(mu_i)
      quantile(mu_i + eps, probs = c(alpha/2, 1 - alpha/2), names = FALSE), numeric(2L)))
    colnames(Q) <- c("Lower_PI","Upper_PI")
    return(Q)
  }
  
  # Draw Theta for Asym / Bimo scenarios
  Theta <- matrix(rnorm(4 * B, mean = 2, sd = 1), ncol = 4, byrow = TRUE)
  MU <- Theta %*% t(X)
  Q <- matrix(NA_real_, nrow = n, ncol = 2)
  
  # Asymmetric gamma errors
  if (scenario == "Asym") {
    EPS <- matrix(rgamma(B * n, shape = 0.1, rate = 0.1), nrow = B, ncol = n)
    
    # Bimodal mixture errors
  } else if (scenario == "Bimo") {
    comp  <- matrix(sample.int(2, size = B * n, replace = TRUE), nrow = B, ncol = n)
    means <- ifelse(comp == 1, -8, 8)
    EPS   <- matrix(rnorm(B * n, mean = as.vector(means), sd = 0.5), nrow = B, ncol = n)
  }
  
  # Compute empirical quantiles
  for (i in 1:n) {
    y <- MU[, i] + EPS[, i]
    Q[i, ] <- quantile(y, probs = c(alpha/2, 1 - alpha/2), names = FALSE)
  }
  colnames(Q) <- c("Lower_PI","Upper_PI")
  Q
}

# ----------------------------------------------
# Run one configuration (n × scenario × missing)
# ----------------------------------------------
run_one_setting <- function(n_val, scenario_val, missing_val,
                            num_simulations = 1000, alpha = 0.1, B_mc = 1e5) {
  n <- n_val; n_test <- 1; m <- rep(5, n); d <- 4
  scenario <- scenario_val
  theta_mu <- rep(2, d); theta_sigma <- diag(1, d)
  X_mu <- rep(0, d); X_sigma <- diag(1, d)
  
  # Missingness mechanism
  if (missing_val == "20") beta <- c(3, 0, 2, 2, 2)
  else beta <- c(0, 0, 2, 2, 2)
  
  # Parallel setup
  num_cores <- parallel::detectCores()
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,
                          c("parent_dir","oracle_pop_interval","B_mc","alpha","scenario"),
                          envir = environment())
  parallel::clusterEvalQ(cl, {
    library(MASS); library(stats); library(lme4); library(merTools)
    source(file.path(parent_dir,"my_functions.R"), chdir = TRUE)
    NULL
  })
  
  # Monte Carlo loop
  results <- foreach::foreach(
    i = 1:num_simulations, .combine = rbind, .options.RNG = 123
  ) %dorng% {
    data <- generate_data(n, m, d, scenario, theta_mu, theta_sigma, X_mu, X_sigma, beta)
    idx <- single_split(n, n_test)
    X_test <- data$X[idx$indices_part1]; Y_test <- data$Y[idx$indices_part1]
    PI <- oracle_pop_interval(as.matrix(X_test[[1]]), scenario, alpha, B_mc, seed = i)
    data.frame(Y_true = Y_test[[1]], Lower_PI_90 = PI[,1], Upper_PI_90 = PI[,2])
  }
  parallel::stopCluster(cl)
  
  # Summary statistics
  covered <- (results$Y_true >= results$Lower_PI_90) &
    (results$Y_true <= results$Upper_PI_90)
  interval_length <- results$Upper_PI_90 - results$Lower_PI_90
  c(mean_coverage = mean(covered), mean_length = mean(interval_length))
}

# ----------------------------------------------
# Main loop over parameter grid
# ----------------------------------------------
param_grid <- expand.grid(
  n = c(100, 300, 500),
  scenario = c("Homo","Heter","Asym","Bimo","Bimo-fix"),
  missing  = c("20","50"),
  stringsAsFactors = FALSE
)

out_list <- vector("list", nrow(param_grid))
for (i in seq_len(nrow(param_grid))) {
  cat(sprintf("\n[%d/%d] n=%s | scenario=%s | missing=%s\n",
              i, nrow(param_grid),
              param_grid$n[i], param_grid$scenario[i], param_grid$missing[i]))
  met <- run_one_setting(param_grid$n[i], param_grid$scenario[i], param_grid$missing[i])
  out_list[[i]] <- data.frame(
    n = param_grid$n[i],
    scenario = param_grid$scenario[i],
    missing = param_grid$missing[i],
    mean_coverage = met["mean_coverage"],
    mean_length   = met["mean_length"]
  )
}

summary_df <- do.call(rbind, out_list)
print(summary_df, row.names = FALSE)

sc_levels <- c("Homo","Heter","Asym","Bimo","Bimo-fix")
summary_df$scenario <- factor(summary_df$scenario, levels = sc_levels)

# Split results into two tables
tab20 <- subset(summary_df, missing == "20",
                select = c(n, scenario, mean_coverage, mean_length))
tab50 <- subset(summary_df, missing == "50",
                select = c(n, scenario, mean_coverage, mean_length))

tab20 <- tab20[order(tab20$n, tab20$scenario), ]
tab50 <- tab50[order(tab50$n, tab50$scenario), ]