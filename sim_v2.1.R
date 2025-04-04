# rm(list = ls())
##########################################################################################
# version 2.1
# data: Apr 4, 2025
# code for data generation
# reference: Yu et al. (2025) Quantile Residual Lifetime Regression for Multivariate Failure Time Data
##########################################################################################

# Load necessary libraries
library(copula)       # For generating dependent errors using Clayton copula
library(survival)     # For Kaplan-Meier estimator
library(quantreg)     # For quantile regression
library(dplyr)        # For data manipulation

# Data generation setting
set.seed(42)
n_sim <- 500
n_cluster <- c(200, 500)
size_cluster <- c(3, 10)
lambda <- 0.69
beta0 <- 1; beta1 <- 1
tau_vals <- c(0.5)
t0_vals <- c(0, 1, 2)
kendalls_taus <- c(0, 0.5, 0.8)
theta_vals <- 2 * kendalls_taus / (1 - kendalls_taus) # Convert Kendall's tau to Clayton theta: theta = 2 * tau / (1 - tau)

all_data <- list() # empty list for data sets
start <- Sys.time()
for(tau in tau_vals){
  for (n in n_cluster) {
    for (m in size_cluster) {
      for (theta in theta_vals) {
        
        kendall_tau <- round(theta / (theta + 2), 2) # convert theta back to Kendall's tau for key setting
        
        scenario_key <- paste0("n", n, "_m", m, "_tau", tau,"_Ktau", kendall_tau)
        
        all_data[[scenario_key]] <- vector(mode = "list", length = n_sim)  # Pre-allocate list
        
        for (b in 1:n_sim) {
          
          x <- matrix(runif(n * m), nrow = n) # individual-level covariate
          
          # copula-based correlated error generation
          # theta <- 8
          # m <- 3
          cop <- claytonCopula(param = theta, dim = m) # theta from theta_vals (line 23) used here
          # p <- rCopula(1000, cop)
          # pairs(p)
          # ?rCopula
          u <- rCopula(n, cop) # generating n rows of m different errors
          eps <- matrix(qexp(u, rate = lambda), nrow = n, ncol = m) 
          
          logT <- beta0 + beta1 * x + log(eps)
          
          T <- exp(logT)
          C <- matrix(runif(n * m, min = 0, max = 20), nrow = n)
          # C <- matrix(data = Inf, nrow = n, ncol = m)
          
          obs_time <- pmin(T, C)
          delta <- (T <= C) * 1
          
          censor_rate <- round(1 - mean(delta), 4)
          
          cat("Sim", b, "n=", n, "m=", m, ": Censoring rate =", censor_rate, "\n")
          
          df <- data.frame(
            id = rep(1:n, each = m),
            x = c(t(x)),
            obs_time = as.vector(obs_time),
            delta = as.vector(delta)
          )
          
          all_data[[scenario_key]][[b]] <- df
          
        }
      }
    }
  }
}
end <- Sys.time()
end - start # 4.97 secs

names(all_data) # synario key list of the generated data set
data <- data.frame(all_data[[1]][1]) # setting
head(data)
1 - mean(data$delta) # censoring rate around 0.3

tau <- 0.5; t0 <- 0
fit <- qris(Surv(obs_time, delta) ~ x, data = data, t0 = t0, Q = tau, nB = 0, method = "nonsmooth")
fit
