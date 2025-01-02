rm(list=ls())
# Load necessary libraries
library(MASS)  # For multivariate normal distribution
library(survival)  # For Kaplan-Meier estimation
library(quantreg)  # For quantile regression

# Data generation function based on cluster-correlated structure
generate_data <- function(cluster_size, num_clusters, rho, censoring_rate, beta) {
  data <- list()
  sigma2 <- 1  # Variance of errors
  for (i in 1:num_clusters) {
    # Covariates
    Z <- rbinom(cluster_size, size = 1, prob = 0.5)
    
    # Errors
    error_cov <- rho * matrix(1, cluster_size, cluster_size) + (1 - rho) * diag(cluster_size)
    e <- MASS::mvrnorm(1, mu = rep(0, cluster_size), Sigma = sigma2 * error_cov)
    
    # Failure times
    lgT <- beta[1] + beta[2] * Z + e  # Linear model
    
    # Censoring times
    C <- rexp(cluster_size, rate = -log(1 - censoring_rate))
    
    # Observed times and indicators
    observed_time <- pmin(lgT, C)
    delta <- as.numeric(lgT <= C)
    
    # Combine data
    cluster_data <- data.frame(cluster = i, Z = Z, T = observed_time, delta = delta)
    data[[i]] <- cluster_data
  }
  return(do.call(rbind, data))
}

# Function to compute L1-type estimation
l1_type_estimation <- function(taus, nsim, cluster_sizes, num_clusters, rho_values, censoring_rates, beta) {
  results <- list()
  for (tau in taus) {
    for (cluster_size in cluster_sizes) {
      for (rho in rho_values) {
        for (censoring_rate in censoring_rates) {
          estimates <- matrix(NA, nrow = nsim, ncol = length(beta))
          for (sim in 1:nsim) {
            # Generate data for the given parameters
            data <- generate_data(cluster_size = cluster_size, num_clusters = num_clusters, rho = rho, censoring_rate = censoring_rate, beta = beta)
            
            # Step 1: Prepare data
            logT <- data$T 
            Z <- data$Z  # Covariates
            
            # Step 2: Kaplan-Meier weights for censoring survival
            km_fit <- survfit(Surv(data$T, 1 - data$delta) ~ 1)
            G_hat <- stepfun(km_fit$time, c(1, km_fit$surv))(data$T)  # Survival function at observed times
            G_hat[G_hat <= 0 | is.na(G_hat)] <- 1  # Handle division by zero or undefined values
            
            # Initialize regression dataset
            U_mat <- cbind(1, Z)
            weights <- data$delta / G_hat  # Inverse probability weights
            
            # Construct pseudo-values for L1 estimation
            pseudo1 <- -apply(U_mat * weights, 2, sum)
            pseudo2 <- 2 * apply(U_mat * weights * tau, 2, sum)
            
            # Combine for the weighted quantile regression system
            Y_reg <- c(logT, rep(0, ncol(U_mat) * 2))
            U_reg <- rbind(U_mat, diag(ncol(U_mat)), diag(ncol(U_mat)))
            wt_reg <- c(weights, rep(1, ncol(U_mat) * 2))
            
            # Fit weighted quantile regression using rq.wfit
            gamma_fit <- tryCatch({
              rq.wfit(U_reg, Y_reg, weights = wt_reg)$coefficients
            }, error = function(e) {
              warning("L1-type estimation failed.")
              return(rep(NA, ncol(U_mat)))
            })
            
            if (!any(is.na(gamma_fit))) {
              estimates[sim, ] <- gamma_fit
            }
          }
          
          # Store averaged results
          key <- paste("Tau:", tau, "Cluster Size:", cluster_size, "Rho:", rho, "Censoring Rate:", censoring_rate)
          results[[key]] <- colMeans(estimates, na.rm = TRUE)
        }
      }
    }
  }
  return(results)
}

# Example Usage
# Set parameter ranges
set.seed(123)
taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
nsim <- 500
cluster_sizes <- c(2, 3) 
num_clusters <- 100 
rho_values <- c(0.2, 0.5, 0.8)
censoring_rates <- c(0.2, 0.4, 0.6) 
beta <- c(2, 1)  # True beta values -- change

# Apply L1-type estimation
start <- Sys.time()
results <- l1_type_estimation(taus, nsim, cluster_sizes, num_clusters, rho_values, censoring_rates, beta)
end <- Sys.time()
end-start

print("Estimated Coefficients for Different Settings:")
print(results)

# # Function to export results to CSV with separate columns for parameters
# export_results_to_csv <- function(results, file_name) {
#   # Initialize an empty data frame to store all results
#   results_df <- data.frame()
#   
#   # Iterate through the list of results
#   for (key in names(results)) {
#     # Extract the corresponding result vector
#     result <- results[[key]]
#     
#     # Parse the key to extract parameter values using regex
#     tau <- as.numeric(sub("Tau: ([0-9.]+).*", "\\1", key))
#     cluster_size <- as.numeric(sub(".*Cluster Size: ([0-9]+).*", "\\1", key))
#     rho <- as.numeric(sub(".*Rho: ([0-9.]+).*", "\\1", key))
#     censoring_rate <- as.numeric(sub(".*Censoring Rate: ([0-9.]+)", "\\1", key))
#     
#     # Create a data frame for the current result
#     temp_df <- data.frame(
#       Tau = tau,
#       Cluster_Size = cluster_size,
#       Rho = rho,
#       Censoring_Rate = censoring_rate,
#       Beta1 = result[1],
#       Beta2 = result[2]
#     )
#     
#     # Append to the results_df
#     results_df <- rbind(results_df, temp_df)
#   }
#   
#   # Export the combined results to a CSV file
#   write.csv(results_df, file = file_name, row.names = FALSE)
# }
# 
# # Use the function to export your results
# export_results_to_csv(results, "estimation_results.csv")

