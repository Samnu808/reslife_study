rm(list = ls())
##########################################################################################
# version 1.2.7
# data: Feb 13, 2025
##########################################################################################
# code for data set generation
##########################################################################################

##########################################################################################
# (I) Load necessary libraries
##########################################################################################
library(copula)
library(MASS)
library(survival)
library(quantreg)

##########################################################################################
# (II) Functions
# (1) gen_copula_err: to generate copula-based errors for clustered survival data
# (2) gen_data: to generate data set
# (3) find_lambda_cens: to find the lambda for censoring time generation 
##########################################################################################

# (II-1) Function to generate copula-based errors for clustered survival data
gen_copula_err <- function(cluster_size, rho, marginal = "normal") {
  
  # Function to convert Kendall's tau (τ) to Pearson's correlation (ρ)
  # The transformation follows the relationship: ρ = sin(π * τ / 2)
  tau_to_rho <- function(tau) sin(pi * tau / 2)
  
  # Convert intra-cluster correlation from Kendall’s tau (τ) to Pearson’s correlation (ρ)
  rho_gaussian <- tau_to_rho(rho)
  
  # Define the copula to generate correlated errors within clusters
  if (marginal == "normal") {
    # Gaussian copula: assumes normally distributed errors with specified correlation
    cop <- normalCopula(rho_gaussian, dim = cluster_size, dispstr = "ex")
  } else if (marginal == "t") {
    # t-Copula: allows for heavier-tailed error distributions (df = 2)
    cop <- tCopula(rho_gaussian, dim = cluster_size, dispstr = "ex", df = 2)
  } else {
    stop("Unsupported marginal distribution! Choose 'normal' or 't'.")
  }
  
  # Generate a single row of correlated uniform random variables from the copula
  u <- rCopula(1, cop)
  
  # Transform uniform samples to the desired marginal distribution
  if (marginal == "normal") {
    errors <- qnorm(u)  # Convert uniform values to standard normal distribution
  } else {
    errors <- qt(u, df = 2)  # Convert uniform values to t-distribution with df = 2
  }
  
  return(as.vector(errors))  # Return the generated vector of correlated errors
}

# (II-2) Function to generate data set
gen_data <- function(cluster_size, num_clusters, tau, rho, beta, marginal, lambda) {
  data_list <- list()  # List to store each cluster's data
  
  # Quantile adjustment (xi: shift parameter)
  xi <- if (tau == 0.5) {
    0
  } else if (marginal == "normal") {
    -qnorm(tau)
  } else {
    -qt(tau, df = 2)
  }
  
  for (i in 1:num_clusters) {  #  Ensure multiple clusters are created
    Z <- rbinom(cluster_size, 1, 0.5)  # Individual-level covariate
    G <- rnorm(1)  # Cluster-level covariate (shared for all failure times in a subject)
    
    # Generate correlated errors for the current cluster
    errors <- gen_copula_err(cluster_size, rho, marginal)
    
    # Add the quantile adjustment
    epsilon <- errors + xi
    
    # Generate log-transformed failure times
    logT <- beta[1] + beta[2] * Z + beta[3] * G + epsilon
    
    # Generate log of censoring times
    if (lambda == 0) {
      logC <- rep(Inf, cluster_size)  # No censoring
    } else {
      logC <- log(rexp(cluster_size, rate = lambda))
    }
    
    # Observed times and event indicators
    log_obs <- pmin(logT, logC)
    delta <- as.numeric(logT <= logC)
    
    # Enforce zero censoring rate when lambda == 0 
    if (lambda == 0) {
      delta <- rep(1, cluster_size)  # Ensure all events are observed
    } else {
      delta[which.max(log_obs)] <- 1  # Ensure the largest observed time is a failure (to make G_hat not 0 later on)
    }
    
    # Store data in a data frame
    cluster_data <- data.frame(
      cluster = rep(i, cluster_size),  #  Multiple clusters now exist
      Z = Z, 
      G = rep(G, cluster_size),
      log_obs = log_obs, 
      delta = delta, 
      logC = logC, 
      logT = logT,
      xi = rep(xi, cluster_size),      # Store xi in each row
      lambda = rep(lambda, cluster_size)  # Store lambda in each row
    )
    
    data_list[[i]] <- cluster_data  #  Each cluster is stored properly
  }
  
  # Combine all clusters into a single dataset
  result <- do.call(rbind, data_list)  #  Ensuring multiple clusters per dataset
  return(result)
}

# (II-3) Function to find the lambda for censoring time generation 
find_lambda_cens <- function(target_rates, taus, rho_values, num_clusters, beta, marginal, max_iter = 100, tolerance = 0.002, initial_lambda = 0.3) {
  results_lambda_rate <- list()
  
  for (tau in taus) {
    for (rho in rho_values) {
      # Create a key for each (tau, rho) combination
      key <- paste0("Tau_", tau, "_Rho_", rho)
      results_lambda_rate[[key]] <- data.frame(lambda = numeric(), target_rate = numeric())  # Keep results as "target_rate"
      
      for (target_rate in target_rates) {
        # Handle target censoring rate = 0
        if (target_rate == 0) {
          lambda <- 0  # No censoring
          data <- gen_data(cluster_size = 2, num_clusters, tau, rho, beta, marginal, lambda)
          achieved_rate <- round(mean(1 - data$delta), 4)  # Achieved censoring rate (should be 0)
          
          # Save results with "target_rate" column name
          results_lambda_rate[[key]] <- rbind(results_lambda_rate[[key]], 
                                              data.frame(lambda = lambda, target_rate = achieved_rate))
          cat("\nTarget Rate is 0: Setting Lambda to 0 | Tau:", tau, "| Rho:", rho, "| Achieved Rate:", achieved_rate, "\n")
          next  # Skip further processing for the target rate = 0
        }
        
        # Default logic for non-zero target rates
        lambda <- initial_lambda  # Initial guess
        prev_lambda <- lambda
        iteration_converged <- TRUE
        oscillation_detected <- FALSE
        
        for (iter in 1:max_iter) {
          # Generate data for a fixed cluster size(=2)
          data <- gen_data(cluster_size = 2, num_clusters, tau, rho, beta, marginal, lambda)
          achieved_rate <- round(mean(1 - data$delta), 4)
          
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
            oscillation_detected <- TRUE
            iteration_converged <- FALSE
            break
          }
        }
        
        # Print the message with "Achieved Rate"
        cat(sprintf("\nTarget Rate is %.1f | Lambda: %.5f | Tau: %.1f | Rho: %.1f | Achieved Rate: %.1f (Iterations: %d)", 
                    target_rate, lambda, tau, rho, achieved_rate, iter), "\n")
        
        # Immediate warnings for issues
        if (oscillation_detected) {
          warning(paste("Lambda oscillation detected for Tau:", tau, ", Rho:", rho, ", Target Rate:", target_rate, ". Stopping early."))
        }
        
        if (!iteration_converged) {
          warning(paste("Iteration did not converge for Tau:", tau, ", Rho:", rho, ", Target Rate:", target_rate))
        }
        
        # Save results for this target rate (keep "target_rate" column)
        results_lambda_rate[[key]] <- rbind(results_lambda_rate[[key]], 
                                            data.frame(lambda = round(lambda, 4), target_rate = achieved_rate))
      }
    }
  }
  
  return(results_lambda_rate)
}

##########################################################################################
# (III) Data generation: parameter settings, etc.
##########################################################################################

# Parameters
set.seed(123)
M <- 1e6
beta <- c(2, 1, 0.5)
taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
nsim <- 500 
cluster_sizes <- c(2, 3)
num_clusters <- 100
rho_values <- c(0.2, 0.5, 0.8)
target_censoring_rates <- c(0, 0.1, 0.3, 0.5)

# Lambda values for the target censoring rates
start <- Sys.time()
lambda_values <- find_lambda_cens(target_censoring_rates, taus, rho_values, num_clusters, beta, "normal")
print(lambda_values)
end <- Sys.time() 
end - start # elapsed time: 1.8 mins

# Extracting and displaying lambda values with all the other parameters
target_censoring_results <- do.call(rbind, lapply(names(lambda_values), function(key) {
  # Extract numeric values for Tau and Rho using improved regex
  tau_value <- as.numeric(gsub("Tau_([0-9.]+)_Rho_.*", "\\1", key))
  rho_value <- as.numeric(gsub(".*_Rho_([0-9.]+)", "\\1", key))
  # Add the Target_Rate column based on predefined target_censoring_rates
  data.frame(
    tau = rep(tau_value, length(lambda_values[[key]]$lambda)),
    rho = rep(rho_value, length(lambda_values[[key]]$lambda)),
    lambda_values[[key]])
}))

# Print the results with the additional Target_Rate column
print(target_censoring_results)

# Generate 500 data sets per setting
data_list <- list()
start <- Sys.time()
for (i in 1:nrow(target_censoring_results)) {
  tau <- target_censoring_results$tau[i]
  rho <- target_censoring_results$rho[i]
  lambda <- target_censoring_results$lambda[i]
  target_rate <- target_censoring_results$target_rate[i]
  
  for (cluster_size in cluster_sizes) {
    # Create a unique key for this setting
    lambda_key <- paste0("Tau_", tau, "_Rho_", rho, "_Cluster_", cluster_size, "_Rate_", target_rate)
    
    # Generate 500 datasets for each (tau, rho, cluster_size, target_rate)
    dataset_list <- vector("list", nsim)  # Pre-allocate list for efficiency
    
    for (sim in 1:nsim) {
      dataset_list[[sim]] <- gen_data(cluster_size, num_clusters, tau, rho, beta, "normal", lambda)
    }
    
    # Store all datasets under the same key
    data_list[[lambda_key]] <- dataset_list
  }
}

end <- Sys.time()
end - start  # Elapsed time: 1.8 hours

# # Check the structure of the generated data_list
# print("Generated data_list keys:")
# print(names(data_list))
# 
# # Preview one dataset to verify correctness
# if (length(data_list) > 0) {
#   print("Preview of the first dataset:")
#   print(head(data_list[[1]]))
# }
# 
# # Check summary statistics of xi and lambda across all datasets
# xi_values <- sapply(data_list, function(df) unique(df$xi))
# lambda_values <- sapply(data_list, function(df) unique(df$lambda))
# 
# print("Unique xi values across datasets:")
# print(unique(xi_values))
# 
# print("Unique lambda values across datasets:")
# print(unique(lambda_values))

# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # Macbook Pro

save(data_list, file = "data_list.RData") # saving the data_list for other simulations

