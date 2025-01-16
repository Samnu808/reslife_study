########### v1.1.1: censoring rates realized are calculated + adjusted parameters are used to achieve the target censoring rates!

rm(list = ls())
# Load necessary libraries
library(copula)    # For copula-based dependency
library(MASS)      # For multivariate normal distribution
library(survival)  # For Kaplan-Meier estimation
library(quantreg)  # For quantile regression
library(ggplot2)   # For plotting boxplots
library(gridExtra) # For arranging multiple plots

# Function to generate covariates
gen_cov <- function(n, k) {
  Z <- matrix(rbinom(n * k, 1, 0.5), nrow = n, ncol = k)
  return(Z)
}

# Function to generate copula-based errors
gen_copula_err <- function(n, k, rho, marginal = "normal") {
  tau_to_rho <- function(tau) sin(pi * tau / 2)
  rho_gaussian <- tau_to_rho(rho)
  
  if (marginal == "normal") {
    cop <- normalCopula(rho_gaussian, dim = k, dispstr = "ex")
  } else if (marginal == "t") {
    cop <- tCopula(rho_gaussian, dim = k, dispstr = "ex", df = 2)
  } else {
    stop("Unsupported marginal distribution! Choose 'normal' or 't'.")
  }
  
  u <- rCopula(n, cop)
  if (marginal == "normal") {
    return(qnorm(u))
  } else if (marginal == "t") {
    return(qt(u, df = 2))
  }
}

# Generate data using copula and xi-based adjustment
gen_data <- function(cluster_size, num_clusters, tau, rho, beta, marginal, censoring_rate) {
  data <- list()
  xi <- if (tau == 0.5) {
    0  # No xi adjustment for tau = 0.5
  } else {
    if (marginal == "normal") {
      -qnorm(tau)
    } else if (marginal == "t") {
      -qt(tau, df = 2)
    }
  }
  
  for (i in 1:num_clusters) {
    Z <- rbinom(cluster_size, 1, 0.5)
    errors <- gen_copula_err(cluster_size, k, rho, marginal)
    epsilon <- errors + xi
    
    # Failure times
    logT <- beta[1] + beta[2] * Z + epsilon
    
    # Censoring times
    if (censoring_rate == 0) {
      C <- rep(Inf, cluster_size)  # No censoring
    } else {
      C <- rexp(cluster_size, rate = -log(1 - censoring_rate))
    }
    
    # Observed times and indicators
    observed_time <- as.numeric(pmin(logT, C))
    delta <- as.numeric(logT <= C)
    
    # Combine data
    cluster_data <- data.frame(cluster = i, Z = Z, T = observed_time, delta = delta, C = C, logT = logT)
    data[[i]] <- cluster_data
  }
  return(do.call(rbind, data))
}

# Adjust censoring rates to match realized censoring rates
adjust_censoring_rate <- function(target_rates, cluster_size, num_clusters, tau, rho, beta, marginal, tolerance = 0.01, max_iter = 100) {
  adjusted_rates <- numeric(length(target_rates))
  adjustment_process <- list()
  
  for (i in seq_along(target_rates)) {
    target_rate <- target_rates[i]
    censoring_rate <- target_rate  # Start with the target rate as the initial guess
    iter <- 1
    diff <- Inf
    iter_log <- data.frame(Iteration = integer(), Censoring_Rate = numeric(), Realized_Rate = numeric())
    
    while (abs(diff) > tolerance && iter <= max_iter) {
      # Generate data with the current censoring_rate
      data <- gen_data(cluster_size, num_clusters, tau, rho, beta, marginal, censoring_rate)
      
      # Calculate the realized censoring rate
      realized_rate <- mean(1 - data$delta)
      diff <- realized_rate - target_rate
      
      # Log the iteration
      iter_log <- rbind(iter_log, data.frame(Iteration = iter, Censoring_Rate = censoring_rate, Realized_Rate = realized_rate))
      
      # Adjust the censoring_rate based on the difference
      censoring_rate <- censoring_rate - diff * 0.1  # Scale factor for adjustment
      censoring_rate <- max(0, min(1, censoring_rate))  # Keep it within [0, 1]
      
      iter <- iter + 1
    }
    
    if (iter > max_iter) {
      warning(paste("Failed to converge for target rate:", target_rate))
    }
    
    adjusted_rates[i] <- censoring_rate
    adjustment_process[[paste0("Target_", target_rate)]] <- iter_log
  }
  
  return(list(Adjusted_Rates = adjusted_rates, Adjustment_Process = adjustment_process))
}

# Function for non-smooth estimation
non_smooth_est <- function(taus, nsim, cluster_sizes, num_clusters, rho_values, censoring_rates, beta, marginal) {
  results <- list()
  beta1_estimates <- data.frame()  # To store beta1 estimates for boxplot
  for (tau in taus) {
    for (cluster_size in cluster_sizes) {
      for (rho in rho_values) {
        for (censoring_rate in censoring_rates) {
          estimates <- matrix(NA, nrow = nsim, ncol = length(beta))
          xi_values <- numeric(nsim)  # Store xi values
          realized_censor_rates <- numeric(nsim)  # Store realized censoring rates
          
          for (sim in 1:nsim) {
            data <- gen_data(cluster_size = cluster_size, num_clusters = num_clusters, tau = tau, 
                             rho = rho, beta = beta, marginal = marginal, censoring_rate = censoring_rate)
            
            logT <- data$T
            Z <- data$Z
            delta <- data$delta
            
            # Calculate realized censoring rate
            realized_censor_rates[sim] <- round(mean(1 - delta), 4)
            
            if (censoring_rate == 0) {
              G_hat <- rep(1, nrow(data))  # No censoring, so weights are all 1
            } else {
              km_fit <- survfit(Surv(data$T, 1 - data$delta) ~ 1)
              G_hat <- stepfun(km_fit$time, c(1, km_fit$surv))(data$T)
              G_hat[G_hat <= 0 | is.na(G_hat)] <- 1
            }
            
            U_mat <- cbind(1, Z)
            weights <- data$delta / G_hat
            
            Y_reg <- c(logT, rep(0, ncol(U_mat) * 2))
            U_reg <- rbind(U_mat, diag(ncol(U_mat)), diag(ncol(U_mat)))
            wt_reg <- c(weights, rep(1, ncol(U_mat) * 2))
            
            gamma_fit <- tryCatch({
              rq.wfit(U_reg, Y_reg, weights = wt_reg)$coefficients
            }, error = function(e) {
              warning("L1-type estimation failed.")
              return(rep(NA, ncol(U_mat)))
            })
            
            if (!any(is.na(gamma_fit))) {
              estimates[sim, ] <- gamma_fit
              xi_values[sim] <- if (tau == 0.5) 0 else {  # Store xi value
                if (marginal == "normal") {
                  -qnorm(tau)
                } else if (marginal == "t") {
                  -qt(tau, df = 2)
                }
              }
              
              # Append beta1 estimate to beta1_estimates
              beta1_estimates <- rbind(beta1_estimates, data.frame(
                Tau = tau,
                Cluster_Size = cluster_size,
                Rho = rho,
                Censoring_Rate = censoring_rate,
                Beta1_Estimate = gamma_fit[2],
                Realized_Censoring_Rate = realized_censor_rates[sim]
              ))
            }
          }
          
          key <- paste("Tau:", tau, "Cluster Size:", cluster_size, "Rho:", rho, "Censoring Rate:", censoring_rate)
          results[[key]] <- list(
            coefficients = colMeans(estimates, na.rm = TRUE),
            xi = mean(xi_values, na.rm = TRUE),
            avg_realized_censor_rate = mean(realized_censor_rates, na.rm = TRUE)
          )
        }
      }
    }
  }
  return(list(results = results, beta1_estimates = beta1_estimates))
}

# Parameters for copula-based data generation and example usage
set.seed(123)
beta <- c(2, 1)  # True beta values
k <- 2            # Number of repeated measures

taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
nsim <- 500
cluster_sizes <- c(2, 3)
num_clusters <- 100
rho_values <- c(0.2, 0.5, 0.8)
target_censoring_rates <- c(0, 0.2, 0.4, 0.6)
marginal <- "normal"  # Change to "t" for t-distribution

# Adjust censoring rates
adjusted_censoring_rates <- adjust_censoring_rate(target_censoring_rates, 2, num_clusters, 0.1, 0.2, beta, marginal)
print(adjusted_censoring_rates)

# Example Usage
start <- Sys.time()
results_list <- non_smooth_est(taus, nsim, cluster_sizes, num_clusters, rho_values, adjusted_censoring_rates$Adjusted_Rates, beta, marginal)
end <- Sys.time()
end - start # elapsed time

results <- results_list$results
beta1_estimates <- results_list$beta1_estimates

# Generate separate boxplots for each censoring rate
gen_plot <- function(tau, cluster_size, rho, data) {
  subset_data <- data[data$Tau == tau & data$Cluster_Size == cluster_size & data$Rho == rho, ]
  plots <- list()
  
  for (censor_rate in unique(subset_data$Censoring_Rate)) {
    censor_subset <- subset_data[subset_data$Censoring_Rate == censor_rate, ]
    
    plot <- ggplot(censor_subset, aes(x = Beta1_Estimate, y = factor(Censoring_Rate))) +
      geom_boxplot() +
      labs(
        x = "Beta1 Estimate",
        y = "Lambda"
      )
    plots[[as.character(censor_rate)]] <- plot
  }
  
  grid.arrange(
    grobs = plots,
    ncol = 1,
    top = paste("Boxplot of Beta1 Estimates for Tau =", tau, ", Cluster Size =", cluster_size, ", Rho =", rho)
  )
}

# Generate boxplots for all combinations
gen_allplots <- function(data) {
  for (tau in unique(data$Tau)) {
    for (cluster_size in unique(data$Cluster_Size)) {
      for (rho in unique(data$Rho)) {
        cat("Generating plots for Tau =", tau, ", Cluster Size =", cluster_size, ", Rho =", rho, "\n")
        gen_plot(tau, cluster_size, rho, data)
      }
    }
  }
}

# Example: Generate all boxplots
gen_plot(0.1,2,0.2,beta1_estimates)
gen_allplots(beta1_estimates)

# Export results to CSV
output <- data.frame()
for (key in names(results)) {
  # Split the key using " " (space) as the delimiter
  tau_cluster_rho_censor <- unlist(strsplit(key, " "))
  
  # Parse values from the split array
  tau <- as.numeric(gsub("Tau:", "", tau_cluster_rho_censor[2]))
  cluster_size <- as.numeric(gsub("Size:", "", tau_cluster_rho_censor[5]))
  rho <- as.numeric(gsub("Rho:", "", tau_cluster_rho_censor[7]))
  censoring_rate <- as.numeric(tau_cluster_rho_censor[10]) # Directly use the numeric value
  
  # Check if values were correctly parsed
  if (is.na(tau) || is.na(cluster_size) || is.na(rho) || is.na(censoring_rate)) {
    warning(paste("Failed to parse key:", key))
    next
  }
  
  # Extract coefficients, xi, and average realized censoring rate
  coef <- results[[key]]$coefficients
  xi <- results[[key]]$xi
  avg_realized_censor_rate <- round(results[[key]]$avg_realized_censor_rate, 4)
  
  # Append to output dataframe
  output <- rbind(output, data.frame(
    Tau = tau,
    Cluster_Size = cluster_size,
    Rho = rho,
    Lambda = censoring_rate,
    Realized_Censor_Rate = avg_realized_censor_rate,
    Beta0 = coef[1],
    Beta1 = coef[2],
    Xi = xi,
    Adjusted_Beta0 = coef[1] - xi
  ))
}

# Write results to CSV
write.csv(output, "results_summary_v1.1.1.csv", row.names = FALSE)
write.csv(beta1_estimates, "beta1_estimates_V1.1.1.csv", row.names = FALSE)



