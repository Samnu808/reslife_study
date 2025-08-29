rm(list = ls()) 
##########################################################################################
# version 2.10
# date: Aug 29, 2025
# code for data generation and simulation: qris function - coef and stderr (is)
#                                          variance estimation w/ partial multiplier bootstrap 
# Scenario 1: individual-level covariate 
# reference: Yu et al. (2025) Quantile Residual Lifetime Regression for Multivariate Failure Time Data
##########################################################################################
# censoring rates: 10% or 30% or 50%
# no coef_bound
# what's added:
# (1) paired mask applied
# (2) cluster-level perturbation: rexp(cluster number)
# (3) boostrap loop: m2 revised
# (4) H, V, Amat: n_cl(number of clusters) instead of N(the total observations) in the denominator
##########################################################################################

# Load required packages
library(survival)
library(quantreg)
library(nleqslv)
library(qris)
library(dplyr)
library(stringr)
library(tidyr)

# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # MacbookPro
load("data_scenario1.RData") # cens30 with maxC=35, beta1=1 (default)
load("data_scenario1_maxC13.RData") # cens50 with maxC=13, beta1=1
load("data_scenario1_maxC90.RData") # cens10 with maxC=90, beta1=1
load("data_scenario1_maxC10.RData") # cens50 with maxC=10, beta1=0
load("data_scenario1_maxC20.RData") # cens30 with maxC=20, beta1=0
load("data_scenario1_maxC60.RData") # cens10 with maxC=60, beta1=0

load("data_scenario2.RData") # cens30 with maxC=35, beta1=1 (default)

# Parameters
beta0 <- 1; beta1 <- 1; p <- 2 # beta1 could be changed
lambda <- 0.69
n_sim <- 100; nb <- 100
M <- 1e6 

# External functions
get_alpha <- function(tau, beta0, beta1, lambda) {
  alpha0 <- log(-log(1 - tau) / lambda) + beta0
  alpha1 <- beta1
  c(alpha0, alpha1)
}

SC.func <- function(T, censor, wgt = 1) {
  deathtime = unique(sort(T[censor == 1]))
  nrisk = ndeath = rep(0, length(deathtime))
  for (i in seq_along(deathtime)) {
    nrisk[i] = sum((deathtime[i] <= T) * wgt)
    ndeath[i] = sum((T == deathtime[i]) * censor * wgt)
  }
  prodobj = 1 - ndeath / nrisk
  survp = sapply(seq_along(prodobj), function(i) prod(prodobj[1:i]))
  return(data.frame(deathtime, ndeath, nrisk, survp))
}

# ---- Paired helpers (mean, sd, mean se, ratio, bias) ----
paired_se_ratio <- function(est_vec, se_vec) {
  mask <- is.finite(est_vec) & is.finite(se_vec)
  if (sum(mask) >= 2L) mean(se_vec[mask]) / sd(est_vec[mask]) else NA_real_
}

paired_stats <- function(est_vec, se_vec, truth_val) {
  mask <- is.finite(est_vec) & is.finite(se_vec)
  if (sum(mask) >= 2L) {
    m_est  <- mean(est_vec[mask])
    sd_est <- sd(est_vec[mask])
    m_se   <- mean(se_vec[mask])
    se_rat <- m_se / sd_est
    bias   <- m_est - truth_val
    c(mean_est = m_est, sd = sd_est, mean_se = m_se,
      se_ratio = se_rat, bias = bias)
  } else {
    c(mean_est = NA_real_, sd = NA_real_, mean_se = NA_real_,
      se_ratio = NA_real_, bias = NA_real_)
  }
}

results <- list()
all_estimates <- list()
all_stderr <- list()

valid_pct_list <- list()
qris_coef_pct_list <- list()
qris_stderr_pct_list <- list()

begin <- Sys.time()

for (scenario_key in names(all_data)) {
  
  tau_val <- as.numeric(str_extract(scenario_key, "(?<=_tau)0\\.\\d+"))
  t0_val <- as.numeric(str_extract(scenario_key, "(?<=_t0)\\d+"))
  truth <- get_alpha(tau_val, beta0, beta1, lambda)
  
  est_matrix <- matrix(NA, nrow = n_sim, ncol = 4)
  stderr_matrix <- matrix(NA, nrow = n_sim, ncol = 4)
  
  colnames(est_matrix)     <- c("IS_0",    "IS_1",    "QRIS_0",    "QRIS_1")
  colnames(stderr_matrix)  <- c("IS_0_se", "IS_1_se", "QRIS_0_se", "QRIS_1_se")
  
  censoring_vec <- numeric(n_sim)
  valid_rate <- numeric(n_sim)
  qris_coef_success_vec <- numeric(n_sim)
  qris_stderr_success_vec <- numeric(n_sim)
  
  for (b in 1:n_sim) {
    
    df <- all_data[[scenario_key]][[b]]
    N <- nrow(df) # from small n to large N (in accordance with Yu et al.)
    n_cl <- length(unique(df$id)) # number of clusters
    censoring_vec[b] <- 1 - mean(df$delta)
    
    # QRIS_IS estimation
    fit_is <- tryCatch({
      qris(Surv(obs_time, delta) ~ x, data = df, t0 = t0_val, Q = tau_val,
           nB = nb, method = "smooth", se = "pmb")
    }, error = function(e) NULL)
    
    if (!is.null(fit_is)) {
      qris_coef_success_vec[b] <- as.numeric(all(is.finite(fit_is$coefficient)))
      qris_stderr_success_vec[b] <- as.numeric(all(is.finite(fit_is$stderr)))
      
      if (qris_coef_success_vec[b] == 1) {
        est_matrix[b, 3:4] <- fit_is$coefficient
      }
      if (qris_stderr_success_vec[b] == 1) {
        stderr_matrix[b, 3:4] <- fit_is$stderr
      }
    }
    
    # Hard-coded IS estimation
    logT <- log(pmax(df$obs_time - t0_val, 1e-4))
    I <- as.numeric(df$obs_time > t0_val)
    # W <- df$delta ### FIX (1): LI'S WEIGHT INSTEAQD OF FAILURE INDICATOR -> W SHOULD BE COMPUTED USIG SC.FUNC (done)
    # Baseline IPCW (unit-weighted KM of censoring)
    
    Gest0 <- SC.func(df$obs_time, 1 - df$delta)
    G_t0 <- if (t0_val > min(Gest0$deathtime)) {
      Gest0$survp[findInterval(t0_val, Gest0$deathtime)]
    } else 1
    idxZ <- findInterval(df$obs_time, Gest0$deathtime)   # may be 0
    G_Z  <- Gest0$survp[idxZ]                            # idx=0 -> NA
    W0   <- df$delta / G_Z * G_t0
    W0[is.na(W0)] <- max(W0, na.rm = TRUE)               # patch NA -> max(W0)
    
    U_mat <- cbind(1, df$x)
    H <- diag(1 / n_cl, p) # N changed to n_cl(number of clusters)
    
    objectF <- function(beta) {
      beta <- as.matrix(beta)
      t(U_mat * I) %*% (W0 * (pnorm((U_mat %*% beta - logT) /
                                      sqrt(diag(U_mat %*% H %*% t(U_mat))))) - tau_val)
    }
    
    # Better initializer (matches qris init="rq")
    beta_init <- tryCatch(
      as.vector(rq.wfit(U_mat, logT, tau = tau_val, weights = I * W0)$coefficients),
      error = function(e) rep(0, p)
    )
    
    is.fit <- tryCatch({
      nleqslv(beta_init, objectF) ### FIX (2): NOT REP(1,P) BUT THE ACTUAL RQ.WFIT RESULTS SHOULD BE INPUT (done)
    }, error = function(e) NULL)
    
    if (!is.null(is.fit) && is.fit$termcd %in% c(1, 2)) {
      est_matrix[b, 1:2] <- is.fit$x
      
      # ISMB Variance estimation
      result.ismb <- c()
      valid_boots <- 0
      
      # outside the bootstrap loop
      m2 <- pnorm((U_mat %*% is.fit$x - logT) / sqrt(diag(U_mat %*% H %*% t(U_mat))))
      
      for (j in 1:nb) {
        
        # cluster-level perturbation 
        id_levels <- sort(unique(df$id))
        id_idx    <- match(df$id, id_levels)      # = as.integer(factor(df$id, levels=id_levels))
        g_i       <- rexp(length(id_levels))      # one gamma per cluster
        eta       <- g_i[id_idx]                  # length N, broadcast per row
        
        # # cluster-level perturbation (not rexp(n))
        # eta <- rexp(n_val)
        # eta <- rep(eta, each = m_val)
        
        W_star <- if (all(df$delta == 1)) {
          rep(1, N)
        } else {
          Gest <- SC.func(df$obs_time, 1 - df$delta, eta)
          # ghat_t0 <- ifelse(t0_val > 0, Gest$survp[findInterval(t0_val, Gest$deathtime)], 1)
          # surv_eta <- Gest$survp[pmin(findInterval(df$obs_time, Gest$deathtime), length(Gest$survp))]
          idx0 <- findInterval(t0_val, Gest$deathtime)
          ghat_t0 <- if (idx0 == 0) 1 else Gest$survp[idx0]
          
          idxZ  <- findInterval(df$obs_time, Gest$deathtime)
          surv_eta <- Gest$survp[idxZ]
          
          W_tmp <- df$delta / surv_eta * ghat_t0
          # W_tmp <- df$delta / Gest$survp[findInterval(df$obs_time, Gest$deathtime)] * ghat_t0
          W_tmp[!is.finite(W_tmp)] <- NA_real_                 # catch Inf/NaN
          W_tmp[is.na(W_tmp)] <- max(W_tmp, na.rm = TRUE)      # NA -> max
          W_tmp
        }
        
        # result <- t(eta * U_mat * I) %*% {
        #   W_star * (pnorm((U_mat %*% is.fit$x - logT) /    # FIX (3): W_star should be put instead of W (done)
        #                sqrt(diag(U_mat %*% H %*% t(U_mat))))) - tau_val
        # }
        
        result <- t(eta * U_mat * I) %*% (W_star * m2 - tau_val) / n_cl
        
        # result <- t(eta * U_mat * I) %*% (
        #   W_star * pnorm((U_mat %*% is.fit$x - logT) /
        #                    sqrt(diag(U_mat %*% H %*% t(U_mat)))) - tau_val
        # ) / N # the same as qris.smooth
        
        if (isTRUE(all(is.finite(result)))) {   # single robust check
          result.ismb <- cbind(result.ismb, result)
          valid_boots <- valid_boots + 1
        }
      }
      
      valid_rate[b] <- valid_boots / nb
      
      v <- cov(t(result.ismb), use = "complete.obs")
      # a.beta <- t(U_mat * I * W0 * as.vector(dnorm((logT - U_mat %*% is.fit$x) /
      #                                               sqrt(diag(U_mat %*% H %*% t(U_mat)))))) %*% (U_mat / sqrt(diag(U_mat %*% H %*% t(U_mat))))
      
      A <- t(U_mat * (I * W0)) %*% (
        - (U_mat / sqrt(diag(U_mat %*% H %*% t(U_mat)))) *
          as.vector(dnorm((U_mat %*% is.fit$x - logT) / sqrt(diag(U_mat %*% H %*% t(U_mat)))))) / n_cl
      A_inv <- qr.solve(A)
      sigma <- t(A_inv) %*% v %*% A_inv
      SE <- sqrt(diag(sigma))
      
      stderr_matrix[b, 1:2] <- SE
    }
    
    # Monitoring output 
    if (b %% n_sim == 0) {
      mean_coef <- colMeans(est_matrix[1:b, ], na.rm = TRUE)
      sd_coef <- apply(est_matrix[1:b, ], 2, sd, na.rm = TRUE)
      mean_se <- colMeans(stderr_matrix[1:b, ], na.rm = TRUE)
      
      # --- REPLACE these two lines with paired versions ---
      # se_ratio <- ifelse(sd_coef < 1e-8, NA, round(mean_se / sd_coef, 2))
      # bias <- round(mean_coef - rep(truth, 2), 4)
      
      # Paired SE ratios
      sr_IS_0   <- paired_se_ratio(est_matrix[1:b, 1], stderr_matrix[1:b, 1])
      sr_IS_1   <- paired_se_ratio(est_matrix[1:b, 2], stderr_matrix[1:b, 2])
      sr_QRIS_0 <- paired_se_ratio(est_matrix[1:b, 3], stderr_matrix[1:b, 3])
      sr_QRIS_1 <- paired_se_ratio(est_matrix[1:b, 4], stderr_matrix[1:b, 4])
      se_ratio <- round(c(sr_IS_0, sr_IS_1, sr_QRIS_0, sr_QRIS_1), 2)
      
      # Paired biases
      pb_IS_0   <- paired_stats(est_matrix[1:b, 1], stderr_matrix[1:b, 1], truth[1])["bias"]
      pb_IS_1   <- paired_stats(est_matrix[1:b, 2], stderr_matrix[1:b, 2], truth[2])["bias"]
      pb_QRIS_0 <- paired_stats(est_matrix[1:b, 3], stderr_matrix[1:b, 3], truth[1])["bias"]
      pb_QRIS_1 <- paired_stats(est_matrix[1:b, 4], stderr_matrix[1:b, 4], truth[2])["bias"]
      bias <- round(c(pb_IS_0, pb_IS_1, pb_QRIS_0, pb_QRIS_1), 4)
      
      mean_cens <- round(mean(censoring_vec[1:b]), 4)
      mean_valid <- round(mean(valid_rate[1:b]), 4)
      
      cat("Sim", b, "|", scenario_key,
          "| Bias:", paste(bias, collapse = " "),
          "| SE_ratio:", paste(se_ratio, collapse = " "),
          "| Censoring:", mean_cens,
          "| Valid Boot %:", mean_valid,
          "| Valid QRIS Coef %:", round(mean(qris_coef_success_vec[1:b], na.rm = TRUE), 4),
          "| Valid QRIS SE %:",   round(mean(qris_stderr_success_vec[1:b], na.rm = TRUE), 4), "\n")
    }
    
  }
  
  # ---- Paired stats per column ----
  ps_IS_0   <- paired_stats(est_matrix[, 1], stderr_matrix[, 1], truth[1])  # hard_IS alpha0
  ps_IS_1   <- paired_stats(est_matrix[, 2], stderr_matrix[, 2], truth[2])  # hard_IS alpha1
  ps_QRIS_0 <- paired_stats(est_matrix[, 3], stderr_matrix[, 3], truth[1])  # qris_smooth alpha0
  ps_QRIS_1 <- paired_stats(est_matrix[, 4], stderr_matrix[, 4], truth[2])  # qris_smooth alpha1
  
  # Save raw matrices (unchanged)
  all_estimates[[scenario_key]] <- as.data.frame(est_matrix)
  all_stderr[[scenario_key]]    <- as.data.frame(stderr_matrix)
  
  # Rates (unchanged)
  valid_pct_list[[scenario_key]]      <- mean(valid_rate, na.rm = TRUE)
  qris_coef_pct_list[[scenario_key]]  <- mean(qris_coef_success_vec, na.rm = TRUE)
  qris_stderr_pct_list[[scenario_key]]<- mean(qris_stderr_success_vec, na.rm = TRUE)
  
  # Assemble summary using paired stats only
  results[[scenario_key]] <- data.frame(
    Coef   = rep(c("alpha0", "alpha1"), 2),
    Method = rep(c("hard_IS", "qris_smooth"), each = 2),
    True   = rep(round(truth, 4), 2),
    
    Est      = round(c(ps_IS_0["mean_est"],  ps_IS_1["mean_est"],
                       ps_QRIS_0["mean_est"],ps_QRIS_1["mean_est"]), 4),
    Bias     = round(c(ps_IS_0["bias"],      ps_IS_1["bias"],
                       ps_QRIS_0["bias"],    ps_QRIS_1["bias"]),     4),
    SD       = round(c(ps_IS_0["sd"],        ps_IS_1["sd"],
                       ps_QRIS_0["sd"],      ps_QRIS_1["sd"]),       4),
    Mean_SE  = round(c(ps_IS_0["mean_se"],   ps_IS_1["mean_se"],
                       ps_QRIS_0["mean_se"], ps_QRIS_1["mean_se"]),  4),
    SE_Ratio = round(c(ps_IS_0["se_ratio"],  ps_IS_1["se_ratio"],
                       ps_QRIS_0["se_ratio"],ps_QRIS_1["se_ratio"]), 4),
    
    Valid_Bootstrap = c(
      rep(round(valid_pct_list[[scenario_key]], 4), 2),
      rep(NA, 2)
    ),
    Valid_QRIS_Coef = c(
      rep(NA, 2),
      rep(round(qris_coef_pct_list[[scenario_key]], 4), 2)
    ),
    Valid_QRIS_SE = c(
      rep(NA, 2),
      rep(round(qris_stderr_pct_list[[scenario_key]], 4), 2)
    )
  )
  
}

end <- Sys.time()

end - begin

warnings()

# Estimation results compiled: results_df
results_df <- bind_rows(results, .id = "Scenario")
print(results_df); str(results_df)

results_df <- results_df %>%
  mutate(
    n   = as.numeric(str_extract(Scenario, "(?<=n)\\d+")),
    m   = as.numeric(str_extract(Scenario, "(?<=_m)\\d+")),
    tau = as.numeric(str_extract(Scenario, "(?<=_tau)\\d*\\.?\\d+")),
    Ktau= as.numeric(str_extract(Scenario, "(?<=_Ktau)\\d*\\.?\\d+")),
    t0  = as.numeric(str_extract(Scenario, "(?<=_t0)\\d+"))
  )

# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # Macbook pro

write.csv(results_df, "SEratio_qrisis_pmb(scenario1,cens30).csv", row.names = FALSE) # Save results_df to CSV
write.csv(results_df, "SEratio_qrisis_pmb(scenario1,cens10,zerobeta1).csv", row.names = FALSE) # Save results_df to CSV

write.csv(results_df, "SEratio_qrisis_pmb(scenario1,cens30,eta_adjusted).csv", row.names = FALSE) # Save results_df to CSV
write.csv(results_df, "SEratio_qrisis_pmb(scenario2,cens30,eta_adjusted).csv", row.names = FALSE) # Save results_df to CSV

est_long <- bind_rows(lapply(names(all_estimates), function(k) {
  df <- all_estimates[[k]]
  df$Scenario <- k
  df$Replication <- 1:nrow(df)
  return(df)
}), .id = "ScenarioID")

est_longer <- est_long %>%
  pivot_longer(
    cols = c("IS_0", "IS_1", "QRIS_0", "QRIS_1"),
    names_to = "Method_Coef",
    values_to = "Estimate"
  ) %>%
  separate(Method_Coef, into = c("Method", "CoefNum"), sep = "_") %>%
  mutate(
    Method = recode(Method, IS = "hard_IS", QRIS = "qris_smooth"),
    Coef = recode(CoefNum, `0` = "alpha0", `1` = "alpha1")
  ) %>%
  mutate(
    n = as.numeric(str_extract(Scenario, "(?<=n)\\d+")),
    m = as.numeric(str_extract(Scenario, "(?<=_m)\\d+")),
    tau = as.numeric(str_extract(Scenario, "(?<=_tau)\\d*\\.?\\d+")),
    Ktau = as.numeric(str_extract(Scenario, "(?<=_Ktau)\\d*\\.?\\d+")),
    t0 = as.numeric(str_extract(Scenario, "(?<=_t0)\\d+"))
  ) %>%
  mutate(True_vals = ifelse(Coef == "alpha0",
                            log(-log(1 - tau) / lambda) + beta0,
                            beta1))

write.csv(est_longer, "SEratio_qrisis_pmb_allest(scenario1,cens30).csv", row.names = FALSE) # Save results_df to CSV
write.csv(est_longer, "SEratio_qrisis_pmb_allest(scenario1,cens10,zerobeta1).csv", row.names = FALSE) # Save results_df to CSV

write.csv(est_longer, "SEratio_qrisis_pmb_allest(scenario1,cens30,eta_adjusted).csv", row.names = FALSE) # Save results_df to CSV
write.csv(est_longer, "SEratio_qrisis_pmb_allest(scenario2,cens30,eta_adjusted).csv", row.names = FALSE) # Save results_df to CSV

se_long <- bind_rows(lapply(names(all_stderr), function(k) {
  df <- all_stderr[[k]]
  df$Scenario <- k
  df$Replication <- 1:nrow(df)
  return(df)
}), .id = "ScenarioID")

write.csv(se_long, "SEratio_qrisis_pmb_allse(scenario1,cens30).csv", row.names = FALSE) # Save results_df to CSV
write.csv(se_long, "SEratio_qrisis_pmb_allse(scenario1,cens10,zerobeta1).csv", row.names = FALSE) # Save results_df to CSV

write.csv(se_long, "SEratio_qrisis_pmb_allse(scenario1,cens30,eta_adjusted).csv", row.names = FALSE) # Save results_df to CSV
write.csv(se_long, "SEratio_qrisis_pmb_allse(scenario2,cens30,eta_adjusted).csv", row.names = FALSE) # Save results_df to CSV

