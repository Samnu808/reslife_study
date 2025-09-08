rm(list = ls()) 
##########################################################################################
# version 2.10
# date: Sep 4, 2025
# code for data generation and simulation: qris function - coef and stderr (is)
#                                          variance estimation w/ partial multiplier bootstrap 
# Scenario 1: individual-level covariate 
# reference: Yu et al. (2025) Quantile Residual Lifetime Regression for Multivariate Failure Time Data
##########################################################################################
# censoring rates: 10% or 30% or 50%
# coef_bound: 10
# what's added: 
# (1) paired mask applied
# (2) cluster-level perturbation: rexp(cluster number)
# (3) boostrap loop: m2 revised
# (4) H, V, Amat: n_cl(number of clusters) instead of N(the total observations) in the denominator
# (5) sel_key introduced: scenarios with only selected parameters can be run
# (6) Plotting code added at the bottom: estimate dist., SE dist., SE_ratio
# (7) Delta(failure indicator) of the largest obs_time manually set to 1 (idx_max delta=1)
# (8) SC.func revision "SC.func_rev": eta ordering changed
# (9) SC.func revision "SC.func_rev2": eta, T, censor ordering changed
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
# load("data_scenario1_maxC13.RData") # cens50 with maxC=13, beta1=1
# load("data_scenario1_maxC90.RData") # cens10 with maxC=90, beta1=1
# load("data_scenario1_maxC10.RData") # cens50 with maxC=10, beta1=0
# load("data_scenario1_maxC20.RData") # cens30 with maxC=20, beta1=0
# load("data_scenario1_maxC60.RData") # cens10 with maxC=60, beta1=0

load("data_scenario2.RData") # cens30 with maxC=35, beta1=1 (default)

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

SC.func_rev <- function(T, censor, wgt = rep(1, length(T))) {
  deathtime = unique(sort(T[censor == 1]))
  wgt = wgt[order(T)] ### added!
  nrisk = ndeath = rep(0, length(deathtime))
  for (i in seq_along(deathtime)) {
    nrisk[i] = sum((deathtime[i] <= T) * wgt)
    ndeath[i] = sum((T == deathtime[i]) * censor * wgt)
  }
  prodobj = 1 - ndeath / nrisk
  survp = sapply(seq_along(prodobj), function(i) prod(prodobj[1:i]))
  return(data.frame(deathtime, ndeath, nrisk, survp))
}

SC.func_rev2 <- function(T, censor, wgt = rep(1, length(T))) {
  idx <- order(T)
  T <- T[idx]; censor <- censor[idx]; wgt <- wgt[idx]   # in accordance with the order of "T"
  deathtime <- unique(T[censor == 1])                   # already sorted in increasing order
  nrisk <- ndeath <- numeric(length(deathtime))
  for (i in seq_along(deathtime)) {
    nrisk[i]  <- sum((deathtime[i] <= T) * wgt)
    ndeath[i] <- sum((T == deathtime[i]) * censor * wgt)
  }
  prodobj <- 1 - ndeath / nrisk
  survp <- cumprod(prodobj)
  data.frame(deathtime, ndeath, nrisk, survp)
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

paired_cov <- function(est, se, truth, z=1.96){
  m <- is.finite(est) & is.finite(se)
  if (sum(m) < 2) return(NA_real_)
  mean(abs(est[m] - truth) <= z * se[m])
}

.kappa_from_eigs <- function(eigs) {
  a <- abs(eigs)
  lam_max <- max(a, na.rm = TRUE)
  lam_min <- max(min(a, na.rm = TRUE), EPS_KAPPA)
  as.numeric(lam_max / lam_min)
}

maybe_print_diag <- function(A, W0, Wstar, tag = "") {
  ## 대칭 Hessian 가정; 대칭이 아닐 수 있으면 symmetric=FALSE로 바꿔도 됨
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  kap <- .kappa_from_eigs(ev)
  if (kap < KAPPA_LOG_THRES) return(invisible(NULL))
  
  qW0  <- unname(quantile(W0,    probs = c(.90, .99, .999), na.rm = TRUE))
  qWst <- unname(quantile(Wstar, probs = c(.90, .99, .999), na.rm = TRUE))
  
  ev_sorted <- sort(ev)
  msg <- sprintf(
    "%s[diag] kappa(A)=%.1f | eig(A)=[%.3g, %.3g] | W0 q90/99/999=%.2f/%.2f/%.2f | W* q90/99/999=%.2f/%.2f/%.2f",
    if (nzchar(tag)) paste0(tag, " ") else "",
    kap, ev_sorted[1], ev_sorted[length(ev_sorted)],
    qW0[1], qW0[2], qW0[3],
    qWst[1], qWst[2], qWst[3]
  )
  cat(msg, "\n")
}

# Parameters
beta0 <- 1; beta1 <- 1 # beta1 could be changed
lambda <- 0.69
n_sim <- 100; nb <- 500
# M <- 1e6 # only for nonsmooth
coef_bound <- 10
t0_vals <- c(0, 1, 2)
taus    <- c(0.25, 0.5, 0.75)

# ---- Floor parameter (single source of truth) ----
EPS_KAPPA <- 1e-12 # kappa 계산용 아주 작은 바닥값(유지)
KAPPA_LOG_THRES <- 40   # kappa(A) 임계치: 이 값 이상일 때만 로그 출력

# Filtering
keep_m    <- c(1, 3, 10)
keep_n    <- c(200, 500)
keep_Ktau <- c(0, 0.5, 0.8)

scenario_keys <- names(all_data)

meta <- tibble::tibble(
  key  = scenario_keys,
  n    = as.numeric(stringr::str_extract(scenario_keys, "(?<=n)\\d+")),
  m    = as.numeric(stringr::str_extract(scenario_keys, "(?<=_m)\\d+")),
  Ktau = as.numeric(stringr::str_extract(scenario_keys, "(?<=_Ktau)\\d+(?:\\.\\d+)?"))
)

sel_keys <- meta %>%
  dplyr::filter(
    n    %in% keep_n,
    m    %in% keep_m,
    Ktau %in% keep_Ktau
  ) %>%
  dplyr::pull(key)

results <- list()
all_estimates <- list()
all_stderr <- list()

valid_pct_list <- list()
qris_coef_pct_list <- list()
qris_stderr_pct_list <- list()

begin <- Sys.time()

for (scenario_key in sel_keys) {
  
  cat(sprintf("\n%s:\n", scenario_key))
  
  for (t0_val in t0_vals) {
    for (tau_val in taus) {
      
      # a unique run key for results
      run_key <- sprintf("%s_tau%.2f_t0%d", scenario_key, tau_val, t0_val)
      
      truth <- get_alpha(tau_val, beta0, beta1, lambda)
      
      est_matrix    <- matrix(NA, nrow = n_sim, ncol = 4)
      stderr_matrix <- matrix(NA, nrow = n_sim, ncol = 4)
      colnames(est_matrix)    <- c("IS_0","IS_1","QRIS_0","QRIS_1")
      colnames(stderr_matrix) <- c("IS_0_se","IS_1_se","QRIS_0_se","QRIS_1_se")
      
      censoring_vec <- numeric(n_sim)
      
      valid_rate <- numeric(n_sim)
      qris_coef_success_vec <- numeric(n_sim)
      qris_stderr_success_vec <- numeric(n_sim)
      
      for (b in 1:n_sim) {
        
        df <- all_data[[scenario_key]][[b]]
        N <- nrow(df) # from small n to large N (in accordance with Yu et al.)
        n_cl <- length(unique(df$id)) # number of clusters
        
        # idx_max <- which.max(df$obs_time) # Identify the index of the maximum obs_time
        # df$delta[idx_max] <- 1 # Force its delta to be 1
        
        # Find all indices with the maximum obs_time
        idx_max_all <- which(df$obs_time == max(df$obs_time, na.rm = TRUE))
        
        # Check for ties
        if (length(idx_max_all) > 1) {
          stop("Multiple maximum obs_time values detected! Simulation stopped.")
        }
        
        # Otherwise, force that one delta to 1
        df$delta[idx_max_all] <- 1L
        
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
        
        # Baseline IPCW via survival::survfit (match QRIS)
        sv <- survival::survfit(Surv(df$obs_time, 1 - df$delta) ~ 1)
        uTime <- sv$time # increasing order (obs_time), showing not just censoring time but all of the obs time
        
        if ( t0_val <= sv$time[1] ) {
          ghat_t0 <- 1
        } else {
          ghat_t0 <- sv$surv[min(which(uTime > t0_val)) - 1] # right before t0_val! (different in bootstrap)
        }
        
        G_Z <- sv$surv[findInterval(df$obs_time, uTime)] # findInterval returns 0 ->
        W0 <- df$delta / G_Z * ghat_t0
        W0[is.na(W0)] <- max(W0, na.rm = TRUE)   # QRIS behavior
        
        # Gest0 <- SC.func(df$obs_time, 1 - df$delta)
        # G_t0 <- if (t0_val > min(Gest0$deathtime)) {
        #   Gest0$survp[findInterval(t0_val, Gest0$deathtime)]
        # } else 1
        # idxZ <- findInterval(df$obs_time, Gest0$deathtime)
        # G_Z_raw <- ifelse(idxZ == 0, 1, Gest0$survp[idxZ]) # S_c(t-) = 1 before first event
        # 
        # floored0_vec[b] <- mean(G_Z_raw <= EPS, na.rm = TRUE)
        # G_Z <- pmax(G_Z_raw, EPS)
        # 
        # W0 <- df$delta / G_Z * G_t0
        # W0[!is.finite(W0)] <- NA_real_ # don't ever try this!!
        # W0[is.na(W0)] <- max(W0, na.rm = TRUE)               # patch NA -> max(W0)
        
        U_mat <- cbind(1, df$x)
        p <- ncol(U_mat)
        
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
        
        if (!is.null(is.fit) && is.fit$termcd %in% c(1, 2) && all(abs(is.fit$x) <= coef_bound)) {
          
          # termcd=1: Function criterion is near zero. Convergence of function values has been achieved.
          # termcd=2: x-values within tolerance. This means that the relative distance between two consecutive x-values is smaller than xtol but that the function value criterion is still larger than ftol. Function values may not be near zero; therefore the user must check if function values are acceptably small.
          est_matrix[b, 1:2] <- is.fit$x
          
          # ISMB Variance estimation
          result.ismb <- NULL
          valid_boots <- 0L
          
          ## inside the replicate, before the bootstrap loop
          alert_count <- 0L
          min_surv_eta <- rep(NA_real_, nb)
          min_idx      <- rep(NA_integer_, nb)  # row index where min happened (optional)
          
          ## before the bootstrap loop
          best_min <- Inf
          
          # outside the bootstrap loop
          m2 <- pnorm((U_mat %*% is.fit$x - logT) / sqrt(diag(U_mat %*% H %*% t(U_mat))))
          
          id_levels <- sort(unique(df$id))
          id_idx    <- match(df$id, id_levels)      # = as.integer(factor(df$id, levels=id_levels))
          
          for (j in 1:nb) {
            
            # cluster-level perturbation 
            g_i       <- rexp(length(id_levels))      # one gamma per cluster
            eta       <- g_i[id_idx]                  # length N, broadcast per row
            
            # # cluster-level perturbation (not rexp(n))
            # eta <- rexp(n_val)
            # eta <- rep(eta, each = m_val)
            
            W_star <- if (all(df$delta == 1)) {
              rep(1, N)
            } else {
              
              sv_w <- survival::survfit(Surv(df$obs_time, 1 - df$delta) ~ 1, weights = eta)
              Gt0_w <- if (t0_val <= sv_w$time[1]) 1 else sv_w$surv[min(which(sv_w$time > t0_val)) - 1]
              idx_w <- findInterval(df$obs_time, sv_w$time)
              GZ_w  <- sv_w$surv[idx_w]
              
              W_tmp <- df$delta / GZ_w * Gt0_w 
              
              # Gest <- SC.func_rev2(df$obs_time, 1 - df$delta, eta)
              # ghat_t0 <- ifelse(t0_val > 0, Gest$survp[findInterval(t0_val, Gest$deathtime)], 1)
              # surv_eta <- Gest$survp[pmin(findInterval(df$obs_time, Gest$deathtime), length(Gest$survp))]
              # idx0 <- findInterval(t0_val, Gest$deathtime)
              # ghat_t0 <- if (idx0 == 0) 1 else Gest$survp[idx0]
              # 
              # idxZ  <- findInterval(df$obs_time, Gest$deathtime)
              # surv_eta <- Gest$survp[idxZ]
              # 
              # W_tmp <- df$delta / surv_eta * ghat_t0
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
            
            if (isTRUE(all(is.finite(result)))) {
              result.ismb <- if (is.null(result.ismb)) result else cbind(result.ismb, result)
              valid_boots <- valid_boots + 1L
            } 
          }
          
          valid_rate[b]   <- valid_boots / nb
          
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
        
        # A의 상태 확인
        condA <- tryCatch(kappa(A), error=function(e) NA_real_)
        eigA  <- tryCatch(range(Re(eigen(A, symmetric=FALSE, only.values=TRUE)$values)), error=function(e) c(NA,NA))
        
        # 가중치 꼬리 확인
        qW0   <- quantile(W0,   c(.90,.99,.999), na.rm=TRUE)
        qWst  <- quantile(W_star,c(.90,.99,.999), na.rm=TRUE)
        
        maybe_print_diag(A, W0, W_star)  # κ(A) ≥ 40 일 때만 한 줄 출력
        
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
          
          mean_cens  <- round(mean(censoring_vec[1:b]), 4)  
          mean_valid <- round(mean(valid_rate[1:b]), 4)
          
          tau_str <- formatC(round(tau_val, 2), format = "f", digits = 2)
          
          cat("Sim", b, "|", paste0("tau=", tau_str, "_t0=", t0_val),
              "| Bias:", paste(bias, collapse=" "),
              "| SE_ratio:", paste(se_ratio, collapse=" "),
              "| Censoring:", mean_cens,
              "| Valid Boot %:",  mean_valid,
              "| Valid QRIS SE %:",   round(mean(qris_stderr_success_vec[1:b], na.rm = TRUE), 4), "\n")
        }
      }
      
      # ---- Paired stats per column ----
      ps_IS_0   <- paired_stats(est_matrix[, 1], stderr_matrix[, 1], truth[1])  # hard_IS alpha0
      ps_IS_1   <- paired_stats(est_matrix[, 2], stderr_matrix[, 2], truth[2])  # hard_IS alpha1
      ps_QRIS_0 <- paired_stats(est_matrix[, 3], stderr_matrix[, 3], truth[1])  # qris_smooth alpha0
      ps_QRIS_1 <- paired_stats(est_matrix[, 4], stderr_matrix[, 4], truth[2])  # qris_smooth alpha1
      
      # Save raw matrices
      all_estimates[[run_key]] <- as.data.frame(est_matrix)
      all_stderr[[run_key]]    <- as.data.frame(stderr_matrix)
      
      # Rates 
      valid_pct_list[[run_key]]      <- mean(valid_rate, na.rm = TRUE)
      qris_coef_pct_list[[run_key]]  <- mean(qris_coef_success_vec, na.rm = TRUE)
      qris_stderr_pct_list[[run_key]]<- mean(qris_stderr_success_vec, na.rm = TRUE)
      
      cov_IS_0   <- paired_cov(est_matrix[,1], stderr_matrix[,1], truth[1])
      cov_IS_1   <- paired_cov(est_matrix[,2], stderr_matrix[,2], truth[2])
      cov_QRIS_0 <- paired_cov(est_matrix[,3], stderr_matrix[,3], truth[1])
      cov_QRIS_1 <- paired_cov(est_matrix[,4], stderr_matrix[,4], truth[2])
      
      # Assemble summary using paired stats only
      results[[run_key]] <- data.frame(
        Coef   = rep(c("alpha0","alpha1"), 2),
        Method = rep(c("hard_IS","qris_smooth"), each = 2),
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
        
        Cov_95   = round(c(cov_IS_0, cov_IS_1, cov_QRIS_0, cov_QRIS_1), 3),
        
        Valid_Bootstrap = c(rep(mean(valid_rate, na.rm=TRUE), 2), rep(NA, 2)),
        Valid_QRIS_Coef = c(rep(NA, 2), rep(mean(qris_coef_success_vec,  na.rm=TRUE), 2)),
        Valid_QRIS_SE   = c(rep(NA, 2), rep(mean(qris_stderr_success_vec, na.rm=TRUE), 2))
      )
    } # end tau loop
  }   # end t0 loop
}     # end scenario_key loop 

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

se_long <- bind_rows(lapply(names(all_stderr), function(k) {
  df <- all_stderr[[k]]
  df$Scenario <- k
  df$Replication <- 1:nrow(df)
  return(df)
}), .id = "ScenarioID")

write.csv(se_long, "SEratio_qrisis_pmb_allse(scenario1,cens30).csv", row.names = FALSE) # Save results_df to CSV
write.csv(se_long, "SEratio_qrisis_pmb_allse(scenario1,cens10,zerobeta1).csv", row.names = FALSE) # Save results_df to CSV

##################################################################################
# Coefficient estimates distribution: boxplots
##################################################################################
# ---- 1) LONG DF FOR ESTIMATES (merge from all_estimates/all_stderr if needed) ----
# all_estimates[[scenario]] has cols: IS_0, IS_1, QRIS_0, QRIS_1
# all_stderr  [[scenario]] has cols: IS_0_se, IS_1_se, QRIS_0_se, QRIS_1_se

est_long <- bind_rows(lapply(names(all_estimates), function(k) {
  df <- all_estimates[[k]]
  df$Scenario <- k
  df$Replication <- seq_len(nrow(df))
  df
}), .id = "ScenarioID")

se_long  <- bind_rows(lapply(names(all_stderr), function(k) {
  df <- all_stderr[[k]]
  df$Scenario <- k
  df$Replication <- seq_len(nrow(df))
  df
}), .id = "ScenarioID")

# Tidy: estimates
est_tidy <- est_long %>%
  pivot_longer(c(IS_0, IS_1, QRIS_0, QRIS_1),
               names_to = "Method_Coef", values_to = "Estimate") %>%
  separate(Method_Coef, into = c("MethodShort","CoefNum"), sep = "_") %>%
  mutate(
    Method = recode(MethodShort, IS = "hard_IS", QRIS = "qris_smooth"),
    Coef   = recode(CoefNum, `0` = "alpha0", `1` = "alpha1")
  ) %>%
  mutate(
    n   = as.numeric(str_extract(Scenario, "(?<=n)\\d+")),
    m   = as.numeric(str_extract(Scenario, "(?<=_m)\\d+")),
    tau = as.numeric(str_extract(Scenario, "(?<=_tau)\\d*\\.?\\d+")),
    Ktau= as.numeric(str_extract(Scenario, "(?<=_Ktau)\\d*\\.?\\d+")),
    t0  = as.numeric(str_extract(Scenario, "(?<=_t0)\\d+"))
  )

# Tidy: standard errors
se_tidy <- se_long %>%
  pivot_longer(c(IS_0_se, IS_1_se, QRIS_0_se, QRIS_1_se),
               names_to = "Method_Coef", values_to = "SE") %>%
  separate(Method_Coef, into = c("MethodShort","CoefNum","se_tag"), sep = "_") %>%
  mutate(
    Method = recode(MethodShort, IS = "hard_IS", QRIS = "qris_smooth"),
    Coef   = recode(CoefNum, `0` = "alpha0", `1` = "alpha1")
  ) %>%
  select(-se_tag)

# Join estimate + SE (useful if you want paired dots/diagnostics later)
est_full <- est_tidy %>%
  left_join(se_tidy, by = c("Scenario","Replication","Method","Coef"))

# Add true values to est_tidy/est_full if you have (beta0, beta1, lambda) in scope
est_full <- est_full %>%
  mutate(True_val = ifelse(Coef == "alpha0",
                           log(-log(1 - tau) / lambda) + beta0,
                           beta1))

# ---- 2) BOXLOT FUNCTION FOR ESTIMATES ----
method_cols <- c("hard_IS" = "#FF6347",  # Tomato
                 "qris_smooth" = "#6495ED") # Cornflower Blue

plot_est_box <- function(data, coef_name, n_val, m_val, tau_val) {
  df <- data %>%
    filter(Coef == coef_name, n == n_val, m == m_val, tau == tau_val,
           is.finite(Estimate))
  
  ggplot(df, aes(x = factor(Ktau), y = Estimate, fill = Method)) +
    geom_boxplot(outlier.size = 0.7) +
    facet_wrap(~ t0, nrow = 1, labeller = labeller(t0 = function(x) paste0("t0 = ", x))) +
    geom_hline(aes(yintercept = True_val), linetype = "dashed", color = "red") +
    labs(
      title = paste0("Estimates of ", coef_name,
                     " | n = ", n_val, ", m = ", m_val,
                     ", tau = ", tau_val),
      x = "Kendall's tau", y = "Estimate"
    ) +
    scale_fill_manual(values = method_cols) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
}

library(ggplot2)

# EXAMPLES:
plot_est_box(est_full, "alpha0", n_val = 200, m_val = 1, tau_val = 0.5)
plot_est_box(est_full, "alpha1", n_val = 200, m_val = 1, tau_val = 0.5)

# Parameters you typically sweep
n_vals   <- sort(unique(est_full$n))
m_vals   <- sort(unique(est_full$m))
tau_vals <- sort(unique(est_full$tau))
coefs    <- c("alpha0", "alpha1")

# Estimates boxplots
for (n0 in n_vals) for (m0 in m_vals) for (tq in tau_vals) for (cc in coefs) {
  p1 <- plot_est_box(est_full, cc, n0, m0, tq)
  ggsave(paste0("box_est_", cc, "_n", n0, "_m", m0, "_tau", tq, ".png"),
         p1, width = 9, height = 5.5, dpi = 300)
  cat("Saved:", paste0("box_est_", cc, "_n", n0, "_m", m0, "_tau", tq, ".png"), "\n")
}

##################################################################################
# Coefficient Standard Error distribution: boxplots
##################################################################################

# all_stderr[[scenario]]: IS_0_se, IS_1_se, QRIS_0_se, QRIS_1_se
se_long <- bind_rows(lapply(names(all_stderr), function(k) {
  df <- all_stderr[[k]]
  df$Scenario <- k
  df$Replication <- seq_len(nrow(df))
  df
}), .id = "ScenarioID")

se_tidy <- se_long %>%
  pivot_longer(c(IS_0_se, IS_1_se, QRIS_0_se, QRIS_1_se),
               names_to = "Method_Coef", values_to = "SE") %>%
  separate(Method_Coef, into = c("MethodShort","CoefNum","se_tag"), sep = "_") %>%
  mutate(
    Method = recode(MethodShort, IS = "hard_IS", QRIS = "qris_smooth"),
    Coef   = recode(CoefNum, `0` = "alpha0", `1` = "alpha1")
  ) %>%
  select(-se_tag) %>%
  mutate(
    n   = as.numeric(str_extract(Scenario, "(?<=n)\\d+")),
    m   = as.numeric(str_extract(Scenario, "(?<=_m)\\d+")),
    tau = as.numeric(str_extract(Scenario, "(?<=_tau)\\d*\\.?\\d+")),
    Ktau= as.numeric(str_extract(Scenario, "(?<=_Ktau)\\d*\\.?\\d+")),
    t0  = as.numeric(str_extract(Scenario, "(?<=_t0)\\d+"))
  ) %>%
  filter(is.finite(SE))

# 색상(이전과 동일 톤)
method_cols <- c("hard_IS" = "orange",   
                 "qris_smooth" = "skyblue")

plot_se_box <- function(data, coef_name, n_val, m_val, tau_val) {
  df <- data %>%
    filter(Coef == coef_name, n == n_val, m == m_val, tau == tau_val)
  
  ggplot(df, aes(x = factor(Ktau), y = SE, fill = Method)) +
    geom_boxplot(outlier.size = 0.7) +
    facet_wrap(~ t0, nrow = 1, labeller = labeller(t0 = function(x) paste0("t0 = ", x))) +
    labs(
      title = paste0("Standard Errors of ", coef_name,
                     " | n = ", n_val, ", m = ", m_val, ", tau = ", tau_val),
      x = "Kendall's tau", y = "Std. Error"
    ) +
    scale_fill_manual(values = method_cols) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
}

plot_se_box(se_tidy, "alpha0", n_val = 200, m_val = 1, tau_val = 0.5)
plot_se_box(se_tidy, "alpha1", n_val = 200, m_val = 1, tau_val = 0.5)

# 배치 저장
n_vals_se   <- sort(unique(se_tidy$n))
m_vals_se   <- sort(unique(se_tidy$m))
tau_vals_se <- sort(unique(se_tidy$tau))
coefs    <- c("alpha0", "alpha1")

for (n0 in n_vals_se) for (m0 in m_vals_se) for (tq in tau_vals_se) for (cc in coefs) {
  p <- plot_se_box(se_tidy, cc, n0, m0, tq)
  ggsave(paste0("box_se_", cc, "_n", n0, "_m", m0, "_tau", tq, ".png"),
         p, width = 9, height = 5.5, dpi = 300)
  cat("Saved:", paste0("box_se_", cc, "_n", n0, "_m", m0, "_tau", tq, ".png"), "\n")
}

plot_se_box_trimmed <- function(data, coef_name, n_val, m_val, tau_val, iqr_mult = 3) {
  df <- data %>%
    filter(Coef == coef_name, n == n_val, m == m_val, tau == tau_val)
  
  # Method × Ktau × t0 그룹별 극단치 제거
  df_trim <- df %>%
    group_by(Method, Ktau, t0) %>%
    mutate(
      Q1 = quantile(SE, 0.01, na.rm = TRUE),
      Q3 = quantile(SE, 0.99, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower = Q1 - iqr_mult * IQR,
      upper = Q3 + iqr_mult * IQR
    ) %>%
    filter(SE >= lower, SE <= upper) %>%
    ungroup()
  
  ggplot(df_trim, aes(x = factor(Ktau), y = SE, fill = Method)) +
    geom_boxplot(outlier.size = 0.7) +
    facet_wrap(~ t0, nrow = 1, labeller = labeller(t0 = function(x) paste0("t0 = ", x))) +
    labs(
      title = paste0("Standard Errors of ", coef_name,
                     " | n = ", n_val, ", m = ", m_val, ", tau = ", tau_val, " (Trimmed)"),
      x = "Kendall's tau", y = "Std. Error"
    ) +
    scale_fill_manual(values = method_cols) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
}

plot_se_box_trimmed(se_tidy, "alpha0", n_val = 200, m_val = 1, tau_val = 0.75)


##################################################################################
# SE ratio: dot plots
##################################################################################
# Assuming results_df has: Tau, n, m, Ktau, t0, Method, SERatio_alpha0, SERatio_alpha1
# Or generically: SERatio, Coef columns
str(results_df)
library(dplyr)
library(ggplot2)

results_df2 <- results_df %>%
  filter(is.finite(SE_Ratio)) %>%
  mutate(
    Method = factor(Method, levels = c("hard_IS", "qris_smooth")),
    Coef   = factor(Coef,   levels = c("alpha0", "alpha1")),
    tau_lab = paste0("tau=", tau),
    t0_lab  = paste0("t0=",  t0)
  )

# First, check the overall SE_Ratio range
range(results_df2$SE_Ratio, na.rm = TRUE)

nm_combos <- results_df2 %>%
  distinct(n, m) %>%
  arrange(n, m)

for (i in 1:nrow(nm_combos)) {
  n_val <- nm_combos$n[i]
  m_val <- nm_combos$m[i]
  
  df_nm <- results_df2 %>% filter(n == n_val, m == m_val)
  
  p <- ggplot(df_nm,
              aes(x = factor(Ktau), y = SE_Ratio,
                  color = Method, shape = Method)) +
    geom_point(size = 2.8, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    facet_grid(rows = vars(Coef), cols = vars(tau_lab, t0_lab)) +
    labs(
      title = paste0("SE Ratio by Ktau | n=", n_val, ", m=", m_val),
      x = "Kendall's tau", y = "SE Ratio"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "right",
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim = c(0.55, 1.25))   # fixed y-axis
  
  # 파일 이름
  fname <- paste0("SEratio_dotplot_n", n_val, "_m", m_val, ".png")
  
  ggsave(fname, plot = p, width = 12, height = 6.5, dpi = 300)
  cat("Saved:", fname, "\n")
}

