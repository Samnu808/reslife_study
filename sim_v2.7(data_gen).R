rm(list = ls())
##########################################################################################
# version 2.7(2.8)
# date: May 22, 2025
# code for data generation 
# Scenario 1: individual-level covariate 
# reference: Yu et al. (2025) Quantile Residual Lifetime Regression for Multivariate Failure Time Data
##########################################################################################
# censoring rates: around 30%

# Load required packages
library(copula)
library(survival)
library(quantreg)
library(emplik)
library(nleqslv)
library(qris)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(forcats)
library(ggh4x)
library(purrr)

##########################################################################################
# STEP 1: Data Generation
##########################################################################################

# Function to compute true alpha coefficients
get_alpha <- function(tau, beta0, beta1, lambda) {
  alpha0 <- log(-log(1 - tau) / lambda) + beta0
  alpha1 <- beta1
  c(alpha0, alpha1)
}

# Parameters 
beta0 <- 1; beta1 <- 1
lambda <- 0.69
M <- 1e6
n_cluster <- c(200, 500)
size_cluster <- c(1, 3, 10)
tau_vals <- c(0.25, 0.5, 0.75)
kendall_taus_all <- c(0, 0.5, 0.8)
t0_vals <- c(0, 1, 2)

set.seed(503)
n_sim <- 500
all_data <- list()

start1 <- Sys.time()

for (n in n_cluster) {
  for (m in size_cluster) {
    for (tau in tau_vals) {
      kendall_taus <- if (m == 1) c(0) else kendall_taus_all
      for (Ktau in kendall_taus) {
        for (t0_gen in t0_vals) {
          key <- paste0("n", n, "_m", m, "_tau", tau, "_Ktau", Ktau, "_t0", t0_gen)
          theta <- 2 * Ktau / (1 - Ktau)
          
          all_data[[key]] <- vector(mode = "list", length = n_sim)
          cumulative_censoring <- numeric(n_sim)
          
          for (b in 1:n_sim) {
            x <- matrix(runif(n * m), nrow = n)
            cop <- claytonCopula(param = theta, dim = m)
            u <- rCopula(n, cop)
            eps_raw <- matrix(qexp(u, rate = lambda), nrow = n, ncol = m)
            
            alpha <- get_alpha(tau, beta0, beta1, lambda)
            logT_raw <- alpha[1] + alpha[2] * x + eps_raw
            T_raw <- t0_gen + exp(logT_raw)
            
            eps_cond <- eps_raw[T_raw > t0_gen]
            q_cond_tau <- quantile(eps_cond, probs = tau)
            eps <- eps_raw - q_cond_tau
            
            logT <- alpha[1] + alpha[2] * x + eps
            T <- t0_gen + exp(logT)
            
            # FIXED maxC = 35 
            C <- matrix(runif(n * m, min = 0, max = 60), nrow = n) # maxC=13(40-50%) or 35(20-30%) or 90(10-20%): beta1=beta0=1
            # maxC=10(40-50%) or 20(20-30%) or 60(10-20%): beta1=0, beta0=1
            obs_time <- pmin(T, C)
            delta <- (T <= C) * 1
            cumulative_censoring[b] <- 1 - mean(delta)
            
            if (b %% n_sim == 0) {
              censor_rate <- round(mean(cumulative_censoring[1:b]), 4)
              cat("Sim", b, key, ": Mean Censoring Rate =", censor_rate, "\n")
            }
            
            df <- data.frame(
              id = rep(1:n, each = m),
              x = as.vector(t(x)),
              obs_time = as.vector(t(obs_time)),
              delta = as.vector(t(delta))
            )
            all_data[[key]][[b]] <- df
          }
        }
      }
    }
  }
}
end1 <- Sys.time()
end1 - start1 # less than 1 min

# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # Macbook Pro

save(all_data, file = "data_scenario1.RData") # saving the data_list for other simulations
save(all_data, file = "data_scenario1_maxC60.RData") # saving the data_list for other simulations
