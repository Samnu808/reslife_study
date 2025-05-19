rm(list = ls())
##########################################################################################
# version 2.6
# date: May 16, 2025
# code for data generation and simulation: non_smooth + induced_smoothing
# Scenario 1: individual-level covariate 
# reference: Yu et al. (2025) Quantile Residual Lifetime Regression for Multivariate Failure Time Data
##########################################################################################
# target censoring rate = 20-30%

# Load required packages
library(copula)
library(survival)
library(quantreg)
library(dplyr)
library(readr)
library(stringr)
library(emplik)
library(nleqslv)
# library(qris)

##########################################################################################
# STEP 1: Data Generation
##########################################################################################

get_alpha <- function(tau, beta0, beta1, lambda) {
  alpha0 <- log(-log(1 - tau) / lambda) + beta0
  alpha1 <- beta1
  c(alpha0, alpha1)
}

# === PARAMETERS ===
beta0 <- 1; beta1 <- 1
lambda <- 0.69
M <- 1e6  # changed from literal value to parameter
n_cluster <- c(200, 500)
size_cluster <- c(1, 3, 10)
tau_vals <- c(0.25, 0.5, 0.75)
kendall_taus_all <- c(0, 0.5, 0.8)
t0_vals <- c(0, 1, 2)
target_rates <- c(0.25)

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
            
            # FIXED maxC = 35 (instead of using maxC_df)
            C <- matrix(runif(n * m, min = 0, max = 35), nrow = n)
            obs_time <- pmin(T, C)
            delta <- (T <= C) * 1
            cumulative_censoring[b] <- 1 - mean(delta)
            
            if (b %% 500 == 0) {
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
end1 - start1

##########################################################################################
# STEP 2: Estimation using fit_non_smooth (qris-style)
##########################################################################################

results <- list()

start2 <- Sys.time()
for (scenario_key in names(all_data)) {
  
  tau_val <- as.numeric(str_extract(scenario_key, "(?<=_tau)0\\.\\d+"))
  t0_val <- as.numeric(str_extract(scenario_key, "(?<=_t0)\\d+"))
  
  est_matrix <- matrix(NA, nrow = n_sim, ncol = 8)  # Now stores both methods
  censoring_vec <- numeric(n_sim)
  
  # Get true alpha values for this scenario
  truth <- get_alpha(tau_val, beta0, beta1, lambda)
  
  for (b in 1:n_sim) {
    df <- all_data[[scenario_key]][[b]]
    
    censoring_vec[b] <- 1 - mean(df$delta)
    
    # qris: non-smooth
    fit_ns <- qris(Surv(obs_time, delta) ~ x, data = df, t0 = t0_val, Q = tau_val, nB = 0, method = "nonsmooth")
    est_matrix[b, 5:6] <- fit_ns$coefficient
    
    fit_is <- qris(Surv(obs_time, delta) ~ x, data = df, t0 = t0_val, Q = tau_val, nB = 0, method = "smooth")
    est_matrix[b, 7:8] <- fit_is$coefficient
    
    logT <- log(pmax(df$obs_time - t0_val, 1e-4)) # Log-transformed observed time, shifted by t0
    I <- as.numeric(df$obs_time > t0_val) # Indicator for being at risk at t0
    
    # Build data frame as qris does
    data <- data.frame(
      Z = df$obs_time,
      logZ = logT,
      I = I,
      delta = df$delta
    )
    
    sv <- survfit(Surv(data$Z, 1 - data$delta) ~ 1) # Kaplan-Meier on all Z
    # ghatt0 as in qris
    if (t0_val <= sv$time[1]) {
      ghatt0 <- 1
    } else {
      ghatt0 <- sv$surv[min(which(sv$time > t0_val)) - 1]
    }
    
    # Weight formula as in qris
    W <- data$delta / sv$surv[findInterval(data$Z, sv$time)] * ghatt0
    W[is.na(W)] <- max(W, na.rm = TRUE)
    data$weight <- W
    
    U_mat <- cbind(1, df$x) # Build U matrix
    
    # Estimating equations using I and W (like qris.nonsmooth)
    pseudo1 <- -colSums(U_mat * I * W)
    pseudo2 <- 2 * colSums(U_mat * I * tau_val)
    
    U_reg <- rbind(U_mat, pseudo1, pseudo2)
    Y_reg <- c(logT, M, M)   # use parameter M instead of 1e6
    wt_reg <- c(W, 1, 1)
    
    fit <- rq.wfit(U_reg, Y_reg, weights = wt_reg)
    est_matrix[b, 1:2] <- fit$coefficients  # Store fit_non_smooth results
    
    # Induced smoothing estimation
    H <- diag(1 / nrow(U_mat), ncol(U_mat), ncol(U_mat))
    
    objectF <- function(beta) {
      beta <- as.matrix(beta)
      result <- t(U_mat * I) %*% (W * (pnorm((U_mat %*% beta - logT) / sqrt(diag(U_mat %*% H %*% t(U_mat)))) - tau_val))
      result <- as.vector(result)
    }
    
    is.fit <- nleqslv(fit$coefficients, objectF)
    if (is.fit$termcd %in% c(1, 2)) {
      est_matrix[b, 3:4] <- is.fit$x  # Store induced smoothing results
    }
    
    if (b %% (n_sim) == 0) {
      mean_fit_ns <- colMeans(est_matrix[1:b, 1:2, drop = FALSE], na.rm = TRUE)
      mean_is <- colMeans(est_matrix[1:b, 3:4, drop = FALSE], na.rm = TRUE)
      mean_ns_qris <- colMeans(est_matrix[1:b, 5:6, drop = FALSE], na.rm = TRUE)
      mean_is_qris <- colMeans(est_matrix[1:b, 7:8, drop = FALSE], na.rm = TRUE)
      
      bias_fit_ns <- round(mean_fit_ns - truth, 4)
      bias_is <- round(mean_is - truth, 4)
      bias_ns_qris <- round(mean_ns_qris - truth, 4)
      bias_is_qris <- round(mean_is_qris - truth, 4)
      
      mean_cens <- round(mean(censoring_vec[1:b]), 4)
      
      cat(
        "Sim", b, "|", scenario_key,
        # "| Mean(NS):", round(mean_fit_ns, 4),
        "| Bias(NS):", bias_fit_ns,
        # "| Mean(IS):", round(mean_is, 4),
        "| Bias(IS):", bias_is,
        # "| Mean(qris_NS):", round(mean_ns_qris, 4),
        "| Bias(qris_NS):", bias_ns_qris,
        # "| Mean(qris_IS):", round(mean_is_qris, 4),
        "| Bias(qris_IS):", bias_is_qris,
        "| Censoring:", mean_cens, "\n")
    }
  }
  
  mean_fit_ns <- colMeans(est_matrix[, 1:2, drop = FALSE], na.rm = TRUE)
  mean_is <- colMeans(est_matrix[, 3:4, drop = FALSE], na.rm = TRUE)
  mean_ns_qris <- colMeans(est_matrix[, 5:6, drop = FALSE], na.rm = TRUE)
  mean_is_qris <- colMeans(est_matrix[, 7:8, drop = FALSE], na.rm = TRUE)
  
  bias_fit_ns <- mean_fit_ns - truth
  bias_is <- mean_is - truth
  bias_ns_qris <- mean_ns_qris - truth
  bias_is_qris <- mean_is_qris - truth
  
  results[[scenario_key]] <- data.frame(
    Coef = rep(c("alpha0", "alpha1"), 4),
    Method = rep(c("fit_non_smooth", "induced_smoothing", "qris_nonsmooth", "qris_smooth"), each = 2),
    True_vals = rep(round(truth,4), 4),
    Est = round(c(mean_fit_ns, mean_is, mean_ns_qris, mean_is_qris), 4),
    Bias = round(c(bias_fit_ns, bias_is, bias_ns_qris, bias_is_qris), 4)
  )
}
end2 <- Sys.time()
end2 - start2 # 500=2.38 hrs, 1000=5.23 hrs

results_df <- bind_rows(results, .id = "Scenario")
print(results_df)

# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # Macbook pro
write.csv(results_df, "results_ns_vs_is_rep1000.csv", row.names = FALSE)


