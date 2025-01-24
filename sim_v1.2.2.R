rm(list=ls())
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
    Z <- rbinom(cluster_size, 1, 0.5)  # Binary covariate
    G <- rnorm(1)  # Standard normal covariate, same for the whole cluster
    errors <- gen_copula_err(cluster_size, 2, rho, marginal)
    epsilon <- errors + xi
    
    # Failure times
    logT <- beta[1] + beta[2] * Z + beta[3] * G + epsilon
    
    # Censoring times
    if (length(censoring_rate) == 1 && censoring_rate == 0) {
      C <- rep(Inf, cluster_size)  # No censoring
    } else {
      C <- rexp(cluster_size, rate = -log(1 - censoring_rate))
    }
    
    # Observed times and indicators
    observed_time <- as.numeric(pmin(logT, C))
    delta <- as.numeric(logT <= C)
    
    # Ensure the largest observed time is not censored
    if (censoring_rate > 0) {
      max_time_idx <- which.max(observed_time)
      delta[max_time_idx] <- 1
    }
    
    # Combine data
    cluster_data <- data.frame(cluster = i, Z = Z, G = rep(G, cluster_size), T = observed_time, delta = delta, C = C)
    data[[i]] <- cluster_data
  }
  return(do.call(rbind, data))
}

# Adjust censoring rates to match realized censoring rates for each tau
adjust_censoring_rate <- function(taus, target_rates, cluster_size, num_clusters, rho, beta, marginal, tolerance = 0.01, max_iter = 100) {
  adjusted_rates <- list()
  adjustment_process <- list()
  
  for (tau in taus) {
    cat("Adjusting censoring rates for tau =", tau, "\n")
    tau_adjusted_rates <- numeric(length(target_rates))
    tau_adjustment_log <- list()
    
    # Calculate xi adjustment for this tau
    xi <- if (tau == 0.5) {
      0  # No xi adjustment for tau = 0.5
    } else {
      if (marginal == "normal") {
        -qnorm(tau)
      } else if (marginal == "t") {
        -qt(tau, df = 2)
      } else {
        stop("Unsupported marginal distribution! Choose 'normal' or 't'.")
      }
    }
    
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
        warning(paste("Failed to converge for target rate:", target_rate, "and tau:", tau))
      }
      
      tau_adjusted_rates[i] <- round(censoring_rate, 4)  # Ensure rounding to 4 decimals
      tau_adjustment_log[[paste0("Target_", target_rate)]] <- iter_log
    }
    
    adjusted_rates[[paste0("Tau_", tau)]] <- tau_adjusted_rates
    adjustment_process[[paste0("Tau_", tau)]] <- tau_adjustment_log
  }
  
  return(list(Adjusted_Rates = adjusted_rates, Adjustment_Process = adjustment_process))
}

# Parameters for copula-based data generation and example usage
set.seed(123)
beta <- c(2, 1, 0.5)  # True beta values
k <- 2                # Number of repeated measures

taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
nsim <- 500
cluster_sizes <- c(2, 3)
num_clusters <- 100
rho_values <- c(0.2, 0.5, 0.8)
target_censoring_rates <- c(0, 0.1, 0.3, 0.5)
marginal <- "normal"  # Change to "t" for t-distribution

# Adjust censoring rates
adjusted_censoring_rates <- adjust_censoring_rate(
  taus = taus,
  target_rates = target_censoring_rates,
  cluster_size = 2,
  num_clusters = num_clusters,
  rho = 0.2,
  beta = beta,
  marginal = marginal
)

print(adjusted_censoring_rates)

# Function for non-smooth estimation
non_smooth_est <- function(taus, nsim, cluster_sizes, num_clusters, rho_values, censoring_rates, beta, marginal) {
  results <- list()
  beta1_estimates <- data.frame()  # To store beta1 estimates for boxplot
  for (tau in taus) {
    for (cluster_size in cluster_sizes) {
      for (rho in rho_values) {
        for (censoring_rate_list in censoring_rates[[paste0("Tau_", tau)]]) {
          estimates <- matrix(NA, nrow = nsim, ncol = length(beta))
          xi_values <- numeric(nsim)  # Store xi values
          realized_censor_rates <- numeric(nsim)  # Store realized censoring rates
          
          for (sim in 1:nsim) {
            data <- gen_data(cluster_size = cluster_size, num_clusters = num_clusters, tau = tau, 
                             rho = rho, beta = beta, marginal = marginal, censoring_rate = censoring_rate_list)
            
            logT <- data$T
            Z <- data$Z
            G <- data$G
            
            # Calculate realized censoring rate
            realized_censor_rates[sim] <- round(mean(1 - data$delta), 4)
            
            if (length(censoring_rate_list) == 1 && censoring_rate_list == 0) {
              G_hat <- rep(1, nrow(data))  # No censoring, so weights are all 1
            } else {
              km_fit <- survfit(Surv(data$T, 1 - data$delta) ~ 1)
              G_hat <- stepfun(km_fit$time, c(1, km_fit$surv))(data$T)
            }
            
            U_mat <- cbind(1, Z, G)
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
                Censoring_Rate = censoring_rate_list,
                Beta1_Estimate = gamma_fit[2],
                Realized_Censor_Rate = realized_censor_rates[sim]
              ))
            }
          }
          
          key <- paste("Tau:", tau, "Cluster Size:", cluster_size, "Rho:", rho, "Censoring Rate:", censoring_rate_list)
          results[[key]] <- list(
            coefficients = round(colMeans(estimates, na.rm = TRUE), 4),
            xi = round(mean(xi_values, na.rm = TRUE), 4),
            avg_realized_censor_rate = round(mean(realized_censor_rates, na.rm = TRUE), 4)
          )
        }
      }
    }
  }
  return(list(results = results, beta1_estimates = beta1_estimates))
}

# Run Simulation
start <- Sys.time()
results_list <- non_smooth_est(taus, nsim, cluster_sizes, num_clusters, rho_values, adjusted_censoring_rates$Adjusted_Rates, beta, marginal)
end <- Sys.time()
end - start # elapsed time

results <- results_list$results
beta1_estimates <- results_list$beta1_estimates

# Function to generate boxplots for Beta1 estimates
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

# Function to generate boxplots for all combinations
gen_allplots <- function(data) {
  for (tau in unique(data$Tau)) {
    for (cluster_size in unique(data$Cluster_Size)) {
      for (rho in unique(data$Rho)) {
        cat("Generating plots for Tau =", tau, ", Cluster Size =", cluster_size, ", Rho =", rho, "\n")
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
    }
  }
}

# Example: Generate all boxplots
gen_plot(0.1, 2, 0.2, beta1_estimates)
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
    Lambda = round(censoring_rate, 4),
    Realized_Censor_Rate = avg_realized_censor_rate,
    Beta0 = coef[1],
    Beta1 = coef[2],
    Beta2 = coef[3],
    Xi = xi,
    Adjusted_Beta0 = coef[1] - xi
  ))
}

setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results")
# Write results to CSV
write.csv(output, "results_summary_v1.2.2.csv", row.names = FALSE)
write.csv(beta1_estimates, "beta1_estimates_v1.2.2.csv", row.names = FALSE)
