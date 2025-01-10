rm(list = ls())
# Load necessary libraries
library(copula)    # For copula-based dependency
library(MASS)      # For multivariate normal distribution
library(survival)  # For Kaplan-Meier estimation
library(quantreg)  # For quantile regression

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
    cluster_data <- data.frame(cluster = i, Z = Z, T = observed_time, delta = delta, C = C)
    data[[i]] <- cluster_data
  }
  return(do.call(rbind, data))
}

# Function for non-smooth estimation
non_smooth_est <- function(taus, nsim, cluster_sizes, num_clusters, rho_values, censoring_rates, beta, marginal) {
  results <- list()
  for (tau in taus) {
    for (cluster_size in cluster_sizes) {
      for (rho in rho_values) {
        for (censoring_rate in censoring_rates) {
          estimates <- matrix(NA, nrow = nsim, ncol = length(beta))
          xi_values <- numeric(nsim)  # Store xi values
          
          for (sim in 1:nsim) {
            data <- gen_data(cluster_size = cluster_size, num_clusters = num_clusters, tau = tau, 
                             rho = rho, beta = beta, marginal = marginal, censoring_rate = censoring_rate)
            
            logT <- data$T
            Z <- data$Z
            
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
            }
          }
          
          key <- paste("Tau:", tau, "Cluster Size:", cluster_size, "Rho:", rho, "Censoring Rate:", censoring_rate)
          results[[key]] <- list(
            coefficients = colMeans(estimates, na.rm = TRUE),
            xi = mean(xi_values, na.rm = TRUE)
          )
        }
      }
    }
  }
  return(results)
}

# Parameters for copula-based data generation and example usage
set.seed(123)
beta <- c(2, 1)  # True beta values
k <- 2            # Number of repeated measures

taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
nsim <- 50
cluster_sizes <- c(2, 3)
num_clusters <- 100
rho_values <- c(0.2, 0.5, 0.8)
censoring_rates <- c(0, 0.2, 0.4, 0.6)
marginal <- "normal"  # Change to "t" for t-distribution

# Example Usage
start <- Sys.time()
results <- non_smooth_est(taus, nsim, cluster_sizes, num_clusters, rho_values, censoring_rates, beta, marginal)
end <- Sys.time() - start

# Export results to CSV
output <- data.frame()
for (key in names(results)) {
  tau_cluster_rho_censor <- unlist(strsplit(key, " "))
  tau <- as.numeric(gsub("Tau:", "", tau_cluster_rho_censor[2]))
  cluster_size <- as.numeric(gsub("Cluster Size:", "", tau_cluster_rho_censor[4]))
  rho <- as.numeric(gsub("Rho:", "", tau_cluster_rho_censor[6]))
  censoring_rate <- as.numeric(gsub("Censoring Rate:", "", tau_cluster_rho_censor[8]))
  
  coef <- results[[key]]$coefficients
  xi <- results[[key]]$xi
  
  output <- rbind(output, data.frame(
    Tau = tau,
    Cluster_Size = cluster_size,
    Rho = rho,
    Censoring_Rate = censoring_rate,
    Beta0 = coef[1],
    Beta1 = coef[2],
    Xi = xi,
    Adjusted_Beta0 = coef[1] - xi
  ))
}

write.csv(output, "results_summary.csv", row.names = FALSE)
