rm(list = ls())
##########################################################################################
# version 2.7.1
# date: May 22, 2025
# code for data generation and simulation: non_smooth + induced_smoothiing
#                                          comparison with qris function
# *code for boxplots added!
# Scenario 2: Cluster-level covariate x_i for all j in cluster i
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
size_cluster <- c(3, 10)
tau_vals <- c(0.25, 0.5, 0.75) # Originally 0.5 only
Ktau <- 0.5  # Fixed for Scenario 2
t0_vals <- c(0, 1, 2)

set.seed(503)
n_sim <- 1000
all_data <- list()

start_gen <- Sys.time()

for (n in n_cluster) {
  for (m in size_cluster) {
    for (tau in tau_vals) {
      for (t0_gen in t0_vals) {
        key <- paste0("n", n, "_m", m, "_tau", tau, "_Ktau", Ktau, "_t0", t0_gen)
        theta <- 2 * Ktau / (1 - Ktau)
        
        all_data[[key]] <- vector(mode = "list", length = n_sim)
        cumulative_censoring <- numeric(n_sim)
        
        for (b in 1:n_sim) {
          x_i <- runif(n)
          x <- matrix(rep(x_i, each = m), nrow = n, ncol = m, byrow = TRUE) # cluster-level covariate
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
          C <- matrix(runif(n * m, min = 0, max = 35), nrow = n)
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
end_gen <- Sys.time()
end_gen - start_gen # 12.6 secs

# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # Macbook Pro

save(all_data, file = "data_scenario2.RData") # saving the data_list for other simulations

##########################################################################################
# STEP 2: Coefficient Estimation (comparing with qris function results)
##########################################################################################
# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # MacbookPro
load("data_scenario2.RData")  # Load pre-generated data set 

results <- list()
all_estimates <- list()
start_sim <- Sys.time()

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
    
    # qris: smooth
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
      result <- t(U_mat * I) %*% (W * (pnorm((U_mat %*% beta - logT) / sqrt(diag(U_mat %*% H %*% t(U_mat))))) - tau_val)
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
        # "| Mean(qris_NS):", round(mean_ns_qris, 4),
        "| Bias(qris_NS):", bias_ns_qris,
        # "| Mean(IS):", round(mean_is, 4),
        "| Bias(IS):", bias_is,
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
  
  all_estimates[[scenario_key]] <- as.data.frame(est_matrix)
  
  results[[scenario_key]] <- data.frame(
    Coef = rep(c("alpha0", "alpha1"), 4),
    Method = rep(c("fit_non_smooth", "qris_nonsmooth", "induced_smoothing", "qris_smooth"), each = 2),
    True_vals = rep(round(truth,4), 4),
    Est = round(c(mean_fit_ns, mean_ns_qris, mean_is, mean_is_qris), 4),
    Bias = round(c(bias_fit_ns, bias_ns_qris, bias_is, bias_is_qris), 4)
  )
}
end_sim <- Sys.time()
end_sim - start_sim # 500=35.2mins

##########################################################################################
# STEP 3: Results - results_df, all_estimates
##########################################################################################

# (1) Estimation results compiled: results_df
results_df <- bind_rows(results, .id = "Scenario")
print(results_df); str(results_df)

# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # Macbook pro
# Save restuls_df to csv file
write.csv(results_df, "results_ns_vs_is_all(scenario2).csv", row.names = FALSE)

# (2) All estimates compiled: all_estimates
est_long <- bind_rows(lapply(names(all_estimates), function(k) {
  df <- all_estimates[[k]]
  colnames(df) <- c("fitNS_0", "fitNS_1", "IS_0", "IS_1", "qrisNS_0", "qrisNS_1", "qrisIS_0", "qrisIS_1")
  df$Scenario <- k
  df$Replication <- 1:nrow(df)
  return(df)
}), .id = "ScenarioID")

str(est_long)

est_longer <- est_long %>%
  pivot_longer(
    cols = all_of(c("fitNS_0", "fitNS_1", "IS_0", "IS_1", 
                    "qrisNS_0", "qrisNS_1", "qrisIS_0", "qrisIS_1")),
    names_to = "Method_Coef",
    values_to = "Estimate"
  ) %>%
  separate(Method_Coef, into = c("Method", "CoefNum"), sep = "_") %>%
  mutate(
    Method = recode(Method, fitNS = "NS", IS = "IS", qrisNS = "qris-NS", qrisIS = "qris-IS"),
    Coef = recode(CoefNum, `0` = "alpha0", `1` = "alpha1")
  )

est_longer <- est_longer %>%
  mutate(
    n = as.numeric(str_extract(Scenario, "(?<=n)\\d+")),
    m = as.numeric(str_extract(Scenario, "(?<=_m)\\d+")),
    tau = as.numeric(str_extract(Scenario, "(?<=_tau)\\d*\\.?\\d+")),
    Ktau = as.numeric(str_extract(Scenario, "(?<=_Ktau)\\d*\\.?\\d+")),
    t0 = as.numeric(str_extract(Scenario, "(?<=_t0)\\d+"))
  )

est_longer <- est_longer %>%
  mutate(True_vals = ifelse(Coef == "alpha0", 
                            log(-log(1 - tau) / lambda) + beta0, 
                            beta1))
# Set Method order
est_longer$Method <- factor(est_longer$Method, levels = c("IS", "qris-IS", "NS", "qris-NS"))

plot_boxp <- function(data, coef_name, cluster_size_val, tau_val, n_val) {
  ggplot(
    data %>% filter(Coef == coef_name, m == cluster_size_val, tau == tau_val, n == n_val),
    aes(x = factor(Ktau), y = Estimate, fill = Method)
  ) +
    geom_boxplot(outlier.size = 0.7) +
    facet_wrap(~ t0, nrow = 1, labeller = labeller(t0 = function(x) paste0("t0 = ", x))) +
    geom_hline(aes(yintercept = True_vals), linetype = "dashed", color = "red") +
    labs(
      title = paste0("Estimates of ", coef_name, 
                     " | n = ", n_val,
                     ", tau = ", tau_val,
                     ", cluster size = ", cluster_size_val),
      x = "Kendall's tau",
      y = "Estimate"
    ) +
    scale_fill_manual(values = c(
      "IS" = "#FF6347",       # Tomato
      "qris-IS" = "#FF8C00",  # Dark Orange
      "NS" = "#6495ED",       # Cornflower Blue
      "qris-NS" = "#9370DB"   # Medium Purple
    )) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

# Draw one plot
plot_boxp(est_longer, coef_name = "alpha0", cluster_size_val = 10, tau_val = 0.25, n_val = 200)
plot_boxp(est_longer, coef_name = "alpha1", cluster_size_val = 3, tau_val = 0.5, n_val = 500)

# Parameters
n_vals <- c(200, 500)
m_vals <- c(3, 10)
tau_vals <- c(0.25, 0.5, 0.75)
coefs <- c("alpha0", "alpha1")

# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # Macbook pro

# Loop to generate and save plots
for (n in n_vals) {
  for (m in m_vals) {
    for (tau in tau_vals) {
      for (coef in coefs) {
        plot_obj <- plot_boxp(est_longer, coef_name = coef, cluster_size_val = m, tau_val = tau, n_val = n)
        
        # Output file name
        file_name <- paste0("boxplot_", coef, "_n", n, "_m", m, "_tau", tau, ".png")
        
        # Save the plot
        ggsave(filename = file_name, plot = plot_obj, width = 9, height = 5.5, dpi = 300)
        
        cat("Saved:", file_name, "\n")
      }
    }
  }
}
