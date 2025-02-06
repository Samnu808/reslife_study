rm(list = ls())

# Load necessary libraries
library(copula)
library(MASS)
library(survival)
library(quantreg)

# Function to generate copula-based errors (corrected to use cluster_size)
gen_copula_err <- function(cluster_size, rho, marginal = "normal") {
  tau_to_rho <- function(tau) sin(pi * tau / 2)
  rho_gaussian <- tau_to_rho(rho)
  
  # Define the copula
  if (marginal == "normal") {
    cop <- normalCopula(rho_gaussian, dim = cluster_size, dispstr = "ex")
  } else if (marginal == "t") {
    cop <- tCopula(rho_gaussian, dim = cluster_size, dispstr = "ex", df = 2)
  } else {
    stop("Unsupported marginal distribution! Choose 'normal' or 't'.")
  }
  
  # Generate one row of correlated errors for the cluster
  u <- rCopula(1, cop)  # Single draw from the copula
  
  if (marginal == "normal") {
    errors <- qnorm(u)  # Transform to normal margins
  } else {
    errors <- qt(u, df = 2)  # Transform to t margins
  }
  
  return(as.vector(errors))  # Return a vector of errors for the cluster
}


# Generate data using copula and xi-based adjustment
gen_data <- function(cluster_size, num_clusters, tau, rho, beta, marginal, lambda) {
  data_list <- list()
  
  # Quantile adjustment
  xi <- if (tau == 0.5) {
    0
  } else if (marginal == "normal") {
    -qnorm(tau)
  } else {
    -qt(tau, df = 2)
  }
  
  for (i in 1:num_clusters) {
    Z <- rbinom(cluster_size, 1, 0.5)  # Failure-time-specific covariate
    G <- rnorm(1)  # Cluster-level covariate (shared for all failure times in subject)
    
    # Generate correlated errors for the current cluster
    errors <- gen_copula_err(cluster_size, rho, marginal)
    
    # Add the quantile adjustment
    epsilon <- errors + xi
    
    # Generate log-transformed failure times
    logT <- beta[1] + beta[2] * Z + beta[3] * G + epsilon
    
    # Generate censoring times
    if (lambda == 0) {
      logC <- rep(Inf, cluster_size)  # No censoring
    } else {
      logC <- log(rexp(cluster_size, rate = lambda))
    }
    
    # Observed times and event indicators
    log_obs <- pmin(logT, logC)
    delta <- as.numeric(logT <= logC)
    delta[which.max(log_obs)] <- 1  # Ensure at least one observed event
    
    # Combine data for the cluster
    cluster_data <- data.frame(
      cluster = rep(i, cluster_size),
      Z = Z, 
      G = rep(G, cluster_size),
      log_obs = log_obs, 
      delta = delta, 
      logC = logC, 
      logT = logT
    )
    data_list[[i]] <- cluster_data
  }
  
  # Combine data from all clusters
  result <- do.call(rbind, data_list)
  return(result)
}

# Function to calculate target censoring rates with adaptive stability measures
target_censoring_rate <- function(target_rates, taus, rho_values, num_clusters, beta, marginal, max_iter = 100, tolerance = 0.002, initial_lambda = 0.3) {
  censoring_results <- list()
  
  for (tau in taus) {
    for (rho in rho_values) {
      # Create a key for each (tau, rho) combination
      key <- paste0("Tau_", tau, "_Rho_", rho)
      censoring_results[[key]] <- data.frame(lambda = numeric(), achieved_rate = numeric())
      
      for (target_rate in target_rates) {
        lambda <- initial_lambda  # Initial guess
        prev_lambda <- lambda
        iteration_converged <- TRUE
        
        for (iter in 1:max_iter) {
          # Generate data for a fixed cluster size (e.g., cluster_size = 2)
          data <- gen_data(cluster_size = 2, num_clusters, tau, rho, beta, marginal, lambda)
          achieved_rate <- round(mean(1 - data$delta), 4)
          
          # Debugging: Print progress
          cat("\nIteration:", iter, "| Tau:", tau, "| Rho:", rho, "| Lambda:", lambda, 
              "| Achieved Rate:", achieved_rate, "| Target Rate:", target_rate, "\n")
          
          # Check stopping criteria
          if (abs(achieved_rate - target_rate) < tolerance) {
            break
          }
          
          # Update lambda proportionally to the error
          adjustment <- (achieved_rate - target_rate) / target_rate
          prev_lambda <- lambda
          lambda <- lambda * (1 - adjustment)
          
          # Dynamically adjust bounds based on the target_rate
          lambda_upper_bound <- max(2, target_rate * 10)  # Empirically determined bound
          lambda <- max(min(lambda, lambda_upper_bound), 1e-4)  # Apply bounds
          
          # Check for oscillation or instability
          if (abs(lambda - prev_lambda) < 1e-6) {
            warning("Lambda oscillation detected. Stopping early.")
            iteration_converged <- FALSE
            break
          }
        }
        
        # Log cases where iteration didn't converge
        if (!iteration_converged) {
          cat("Warning: Iteration did not converge for Tau:", tau, ", Rho:", rho, 
              ", Target Rate:", target_rate, "\n")
        }
        
        # Save results for this target rate
        censoring_results[[key]] <- rbind(censoring_results[[key]], 
                                          data.frame(lambda = round(lambda, 4), achieved_rate = achieved_rate))
      }
    }
  }
  
  return(censoring_results)
}

# Non-smooth estimation function
non_smooth_est <- function(taus, nsim, cluster_sizes, num_clusters, rho_values, censoring_results, beta, marginal) {
  results <- list()
  
  for (tau in taus) {
    for (rho in rho_values) {
      for (cluster_size in cluster_sizes) {
        # Construct the key for accessing censoring rates
        censoring_lambda <- paste0("Tau_", tau, "_Rho_", rho)
        
        # Check if the key exists in the censoring_results list
        if (!censoring_lambda %in% names(censoring_results)) {
          warning("Censoring rates key not found:", censoring_lambda)
          next
        }
        
        # Retrieve the lambda values for the given tau and rho
        lambda_values <- censoring_results[[censoring_lambda]]$lambda
        
        # Loop through the lambda values for the current (tau, rho) configuration
        for (lambda in lambda_values) {
          estimates <- matrix(NA, nrow = nsim, ncol = length(beta))
          xi_values <- numeric(nsim)
          realized_censor_rates <- numeric(nsim)
          
          for (sim in 1:nsim) {
            # Generate data for the current simulation
            data <- gen_data(cluster_size, num_clusters, tau, rho, beta, marginal, lambda)
            logT <- data$log_obs  
            delta <- data$delta
            Z <- data$Z
            G <- data$G
            
            # Calculate realized censoring rate
            realized_censor_rates[sim] <- mean(1 - delta)
            
            # Calculate G_hat for survival probabilities
            if (all(delta == 1)) {
              G_hat <- rep(1, length(logT))
            } else {
              km_fit <- survfit(Surv(logT, 1 - delta) ~ 1)
              G_hat <- stepfun(km_fit$time, c(1, km_fit$surv))(logT)
            }
            
            # Construct the design matrix and calculate weights
            U_mat <- cbind(1, Z, G)
            weights <- delta / pmax(G_hat, 1e-6)  # Avoid division by zero
            
            # Calculate pseudo-responses for L1-type estimation
            pseudo1 <- -apply(U_mat * weights, 2, sum)
            pseudo2 <- 2 * apply(U_mat * weights * tau, 2, sum)
            
            # Construct the regression matrix and response vector
            U_reg <- rbind(U_mat, pseudo1, pseudo2)
            Y_reg <- c(logT, M, M)
            wt_reg <- c(weights, 1, 1)
            
            # Perform quantile regression using L1 minimization
            gamma_fit <- tryCatch({
              rq.wfit(U_reg, Y_reg, weights = wt_reg)$coefficients
            }, error = function(e) {
              warning("rq.wfit failed:", conditionMessage(e))
              return(rep(NA, ncol(U_mat)))
            })
            
            # Store the results if regression was successful
            if (!any(is.na(gamma_fit))) {
              estimates[sim, ] <- gamma_fit
              if (tau == 0.5) {
                xi_values[sim] <- 0
              } else if (marginal == "normal") {
                xi_values[sim] <- -qnorm(tau)
              } else {
                xi_values[sim] <- -qt(tau, df = 2)
              }
            }
          }
          
          # Store results for the current parameter combination
          key <- paste("Tau:", tau, "Rho:", rho, "Cluster Size:", cluster_size, "Lambda:", lambda)
          results[[key]] <- list(
            coefficients = colMeans(estimates, na.rm = TRUE),
            xi = mean(xi_values, na.rm = TRUE),
            avg_realized_censor_rate = mean(realized_censor_rates, na.rm = TRUE)
          )
        }
      }
    }
  }
  
  return(results)
}

# Simulation parameters
set.seed(123)
beta <- c(2, 1, 0.5)
M <- 1e6 

taus <- c(0.1,0.3,0.5,0.9)
# taus <- c(0.1)
nsim <- 500
cluster_sizes <- c(2, 3)
num_clusters <- 100
rho_values <- c(0.2,0.5,0.8)
# rho_values <- c(0.2)
target_censoring_rates <- c(0, 0.1, 0.3, 0.5)

# lambda values
target_censoring_results <- target_censoring_rate(target_censoring_rates, taus, rho_values, num_clusters, beta, "normal")
print(target_censoring_results)

# Run the simulation
start_time <- Sys.time()
results_list <- non_smooth_est(taus, nsim, cluster_sizes, num_clusters, rho_values, target_censoring_results, beta, "normal")
end_time <- Sys.time()
end_time - start_time

# warnings()

############################################################################################################################

# Initialize an empty data frame to store the results
output <- data.frame()

# Iterate over each key in the results_list
for (key in names(results_list)) {
  # Split the key to extract parameters
  params <- unlist(strsplit(key, " "))
  tau <- as.numeric(params[2])
  rho <- as.numeric(params[4])
  cluster_size <- as.numeric(params[7])
  lambda <- as.numeric(gsub(",", ".", params[9]))  # Handle numeric conversion
  
  # Retrieve the estimated coefficients and other metrics
  coef <- results_list[[key]]$coefficients
  xi <- results_list[[key]]$xi
  avg_realized_censor_rate <- results_list[[key]]$avg_realized_censor_rate
  
  # Append the extracted and rounded information to the output data frame
  output <- rbind(output, data.frame(
    Tau = tau,
    Cluster_Size = cluster_size,
    Rho = rho,
    Lambda = format(round(lambda, 4), nsmall = 4, scientific = FALSE),
    Censoring_Rate = format(round(avg_realized_censor_rate, 4), nsmall = 4, scientific = FALSE),
    Beta0 = format(round(coef[1], 4), nsmall = 4, scientific = FALSE),
    Beta1 = format(round(coef[2], 4), nsmall = 4, scientific = FALSE),
    Beta2 = format(round(coef[3], 4), nsmall = 4, scientific = FALSE),
    Xi = format(round(xi, 4), nsmall = 4)
  ))
}

# Check the structure of the output data frame
print(output)

# Write results to CSV
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results")
write.csv(output, "results_summary_v1.2.4.csv", row.names = FALSE)
