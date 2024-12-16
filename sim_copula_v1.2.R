
#################################################################################
### Code for simulation: data generation by copulas (Fu and Wang, 2016)
#################################################################################
### Version 1.2
### Date: Dec-17-2024 
#################################################################################
### Changes from previous version: 
### (1) more various tau's and rho's 
### (2) adoption of xi_initial and xi_adjusted 
### (3) Cases 1(Multivariate normal) and 4(Non-normal: multivariate-t) 
#################################################################################
### Contents: 
### (I) Required packages and key parameter settings
### (II) Functions for data generation 
### (III) Simulation: Cases 1 and 4 
### (IV) Result processing and documentation 
#################################################################################

rm(list = ls())
#################################################################################
### (I) Required packages and key parameter settings
#################################################################################
# Load required libraries
library(copula)    # For copula-based dependency
library(quantreg)  # For quantile regression

# Parameters 
set.seed(1)
n <- 150          # Number of subjects
k <- 6            # Number of repeated measures
beta <- c(1, 1)   # True coefficients
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)  # Quantiles
rho_vals <- c(0, 0.2, 0.5, 0.8)       # Kendall's tau for copula dependency
reps <- 1000      # Number of simulation repetitions


#################################################################################
### (II) Functions for data generation 
#################################################################################
# Generate covariates
gen_cov <- function(n, k) {
  x1 <- matrix(rbinom(n * k, 1, 0.5), nrow = n, ncol = k) # Binary covariate
  x2 <- matrix(rep(1:k, each = n), nrow = n, ncol = k)    # Observation times
  list(x1 = x1, x2 = x2)
}

# Generate copula-based multivariate errors
gen_copula_err <- function(n, k, rho, marginal = "normal") {
  # Convert Kendall's tau to Pearson's correlation for Gaussian copula
  tau_to_rho <- function(tau) sin(pi * tau / 2)
  rho_gaussian <- tau_to_rho(rho)
  
  if (marginal == "normal") {
    cop <- normalCopula(rho_gaussian, dim = k, dispstr = "ex")
  } else if (marginal == "t") {
    cop <- tCopula(rho_gaussian, dim = k, dispstr = "ex", df = 2)
  } else {
    stop("Unsupported marginal distribution! Choose 'normal' or 't'.")
  }
  
  # Generate copula-based uniform random variables
  u <- rCopula(n, cop)
  
  # Transform to desired marginal distributions
  if (marginal == "normal") {
    qnorm(u, mean = 0, sd = 1)
  } else if (marginal == "t") {
    qt(u, df = 2)
  }
}

# Generate data for Cases 1(MVN) and 4(multivariate-t)
gen_data <- function(n, k, beta, x1, x2, tau, rho, xi, marginal) {
  errors <- gen_copula_err(n, k, rho, marginal)
  epsilon <- errors + xi
  y <- beta[1] * x1 + beta[2] * x2 + epsilon
  list(y = y, x1 = x1, x2 = x2)
}


#################################################################################
### (III) Simulation: Cases 1 and 4 
#################################################################################
results_copula <- list() # Store the simulation results
begin <- Sys.time() # Record the current time at the beginning

for (marginal in c("normal", "t")) {
  for (tau in taus) {
    for (rho in rho_vals) {
      covariates <- gen_cov(n, k)
      
      if(marginal == "normal"){
        xi_initial <- -qnorm(tau) # Initial guess for xi (Case 1)
      } else if (marginal == "t"){
        xi_initial <- -qt(tau, df=2) # Initial guess for xi (Case 4)
      }
      
      # Perform all repetitions without iterative adjustment
      sim_results <- replicate(reps, {
        data <- gen_data(n, k, beta, covariates$x1, covariates$x2, tau, rho, xi_initial, marginal)
        y_vec <- as.vector(data$y)
        x1_vec <- as.vector(data$x1)
        x2_vec <- as.vector(data$x2)
        X <- cbind(x1_vec, x2_vec)
        fit <- rq(y_vec ~ X, tau = tau)
        coef(fit)
      })
      
      # Compute adjusted xi
      intercept_mean <- mean(sim_results[1, ])  # Mean of the intercept term (should be zero)
      xi_adjusted <- xi_initial - intercept_mean 
      
      # Store results
      results_copula[[paste0("tau=", tau, "_rho=", rho, "_marginal=", marginal)]] <- list(
        summary = apply(sim_results, 1, function(x) c(mean = mean(x), sd = sd(x))),
        xi_initial = xi_initial,
        xi_adjusted = xi_adjusted
      )
    }
  }
}

end <- Sys.time() # Record the time at the end of the simulation

elapsed_t <- end - begin; print(elapsed_t) # Check the duration of the process

#################################################################################
### (IV) Result processing and documentation 
#################################################################################
# Produce a summary table based on the simulation 
output_table <- do.call(rbind, lapply(names(results_copula), function(name) {
  res <- results_copula[[name]]
  summary <- res$summary
  
  # Debugging: Check structure of summary
  print(paste("Case:", name))
  print(summary)
  
  # Correct extraction of row/column data based on structure of `summary`
  data.frame(
    Case = name,
    xi_initial = round(res$xi_initial,6),
    xi_adjusted = round(res$xi_adjusted,6),
    Intercept_mean = round(summary["mean", "(Intercept)"],6),
    Intercept_sd = round(summary["sd", "(Intercept)"],6),
    x1_mean = round(summary["mean", "Xx1_vec"],6),
    x1_sd = round(summary["sd", "Xx1_vec"],6),
    x2_mean = round(summary["mean", "Xx2_vec"],6),
    x2_sd = round(summary["sd", "Xx2_vec"],6)
  )
}))

print(output_table) # Print Results

# Save the table to CSV
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results")
write.csv(output_table, paste0("results_",as.character(Sys.Date()),".csv"), row.names = FALSE)

