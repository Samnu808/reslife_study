rm(list = ls())
##########################################################################################
# version 1.2.9.2
# data: Mar 27, 2025
# modifications from earlier version:
# (1) compare_methods function: warning message capturing part is revised 
# (2) boxp_combined function is added for diff data sets comparison
##########################################################################################

##########################################################################################
# (I) Load necessary libraries
##########################################################################################
# packages to load
library(survival)
library(quantreg)
library(qris)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyr)
library(crayon) # mutate

##########################################################################################
# (II) Functions
# (1) capt_warnin: Function to capture warnings and store detailed information
# (2) compare_methods: to compare two methods(qris vs. crq vs. Non-Smooth Estimation)
##########################################################################################
# (II-1) Function to capture warnings and store detailed information
warning_log <- list() # List to store warnings
capt_warnin <- function(tau, rho, cluster_size, rate, method, message) {
  warning_log <<- append(warning_log, list(
    list(Tau = tau, Rho = rho, Cluster_Size = cluster_size, Target_Rate = rate, Method = method, Message = message)
  ))
}
# (II-2) Function to compare two methods: QRIS and Non-Smooth Estimation
compare_methods <- function(data_list, nsim, beta, taus, rho_values, cluster_sizes, target_censoring_rates) {
  results_df <- data.frame()
  beta_estimates_df <- data.frame()
  
  warning_log <<- data.frame(Tau = numeric(), Rho = numeric(), Cluster_Size = numeric(), 
                             Target_Rate = numeric(), Method = character(), Message = character(),
                             stringsAsFactors = FALSE)
  
  for (tau in taus) {
    for (rho in rho_values) {
      for (cluster_size in cluster_sizes) {
        for (rate in target_censoring_rates) {
          
          rate_str <- format(round(rate, 1), nsmall = 1)
          lambda_key <- paste0("Tau_", tau, "_Rho_", rho, "_Cluster_", cluster_size, "_Rate_", rate_str)
          
          if (!lambda_key %in% names(data_list)) next
          
          datasets <- data_list[[lambda_key]]  
          realized_rates <- numeric(nsim)
          
          convergence_count_qris <- 0
          convergence_count_crq <- 0
          convergence_count_non_smooth <- 0
          nonunique_count_qris <- 0
          nonunique_count_crq <- 0
          nonunique_count_non_smooth <- 0  
          
          beta_qris <- matrix(NA, nrow = nsim, ncol = length(beta))
          beta_crq <- matrix(NA, nrow = nsim, ncol = length(beta))
          beta_non_smooth <- matrix(NA, nrow = nsim, ncol = length(beta))
          
          xi_val <- unique(datasets[[1]]$xi)[1]
          
          for (sim in 1:nsim) {
            if (sim %% msg_freq == 0) {
              cat("\n--- Iteration:", sim, "of", nsim, "completed for", lambda_key, "---\n")
            }
            
            data <- datasets[[sim]]
            data <- data[order(data$log_obs), ]
            lambda_value <- unique(data$lambda)
            
            logT <- data$log_obs
            delta <- data$delta
            Z <- data$Z
            G <- data$G
            
            realized_rates[sim] <- mean(1 - delta)
            
            sv <- survfit(Surv(logT, 1 - delta) ~ 1)
            t0 <- 0
            
            if (t0 <= sv$time[1]) {
              ghatt0 <- 1
            } else {
              ghatt0 <- sv$surv[min(which(sv$time > t0)) - 1]
            }
            
            if (all(delta == 1)) {
              # G_hat <- rep(1, length(logT)) 
              W <- rep(1, length(logT))  # Use same logic as QRIS
              
            } else {
              # km_fit <- survfit(Surv(logT, 1 - delta) ~ 1)
              # G_hat <- stepfun(km_fit$time, c(1, km_fit$surv))(logT)
              
              W <- delta / sv$surv[findInterval(logT, sv$time)] * ghatt0
              W[is.na(W)] <- max(W, na.rm = TRUE)
            }
            
            U_mat <- cbind(1, Z, G)
            
            pseudo1 <- -apply(U_mat * W, 2, sum)
            pseudo2 <- 2 * apply(U_mat * tau, 2, sum)
            
            U_reg <- rbind(U_mat, pseudo1, pseudo2)
            Y_reg <- c(logT, M, M)
            wt_reg <- c(W, 1, 1)
            
            # qris model fitting
            fit_qris <- tryCatch(
              withCallingHandlers(
                {
                  qris(Surv(exp(log_obs), delta) ~ Z + G, data = data, t0 = 0, Q = tau, nB = 0, method = "nonsmooth")
                },
                warning = function(w) {
                  capt_warnin(tau, rho, cluster_size, rate, "QRIS", w$message)
                  if (grepl("nonunique", w$message, ignore.case = TRUE)) {
                    nonunique_count_qris <<- nonunique_count_qris + 1
                  }
                  invokeRestart("muffleWarning")
                }
              ),
              error = function(e) {
                capt_warnin(tau, rho, cluster_size, rate, "QRIS", e$message)
                return(NULL)
              }
            )
            
            # crq fitting
            tau_grid <- c(tau - 0.1, tau, tau + 0.1) # used only for crq function
            fit_crq <- tryCatch(
              {
                fit <- crq(Surv(log_obs, delta) ~ Z + G, data = data, method = "Portnoy")
                summary_fit <- summary(fit, taus = tau_grid)
                
                if (length(summary_fit) > 0) {
                  coef_crq <- summary_fit[[which(round(sapply(summary_fit, function(x) x$tau), 1) == tau)]]$coefficients[, "Value"]
                  convergence_count_crq <- convergence_count_crq + 1
                  coef_crq  
                } else {
                  rep(NA, length(beta))
                }
              },
              warning = function(w) {
                capt_warnin(tau, rho, cluster_size, rate, "CRQ", w$message)
                nonunique_count_crq <- nonunique_count_crq + 1  
                return(rep(NA, length(beta)))
              },
              error = function(e) return(rep(NA, length(beta)))
            )
            
            # Non-Smooth Estimation 
            fit_non_smooth <- tryCatch(
              withCallingHandlers(
                {
                  result <- rq.wfit(U_reg, Y_reg, weights = wt_reg)$coefficients
                  if (any(is.na(result))) {
                    warning("Non-smooth fit resulted in NA coefficients.")
                  }
                  result
                },
                warning = function(w) {
                  capt_warnin(tau, rho, cluster_size, rate, "NonSmooth", w$message)
                  if (grepl("nonunique", w$message, ignore.case = TRUE)) {
                    nonunique_count_non_smooth <<- nonunique_count_non_smooth + 1
                  }
                  invokeRestart("muffleWarning")
                }
              ),
              error = function(e) {
                capt_warnin(tau, rho, cluster_size, rate, "NonSmooth", e$message)
                return(rep(NA, length(beta)))
              }
            )
            
            # Store Results Only if Fit Exists
            if (!is.null(fit_qris)) {
              beta_qris[sim, ] <- coef(fit_qris)  # No unnecessary modifications to coef()
              convergence_count_qris <- convergence_count_qris + 1
            } else {
              beta_qris[sim, ] <- rep(NA, length(beta))
            }
            
            if (!is.null(fit_crq) && !all(is.na(fit_crq))) {
              beta_crq[sim, ] <- fit_crq
              convergence_count_crq <- convergence_count_crq + 1
            } else {
              beta_crq[sim, ] <- rep(NA, length(beta))
            }
            
            if (!is.null(fit_non_smooth) && !all(is.na(fit_non_smooth))) {
              beta_non_smooth[sim, ] <- fit_non_smooth
              convergence_count_non_smooth <- convergence_count_non_smooth + 1
            } else {
              beta_non_smooth[sim, ] <- rep(NA, length(beta))
            }
            
            if (sim %% msg_freq == 0) {
              cat("Cumulative Mean Beta (QRIS) at", sim, "iterations:","\n")
              print(round(colMeans(beta_qris[1:sim, ], na.rm = TRUE), 4))
              
              cat("Cumulative Mean Beta (CRQ) at", sim, "iterations:","\n")
              print(round(colMeans(beta_crq[1:sim, ], na.rm = TRUE), 4))
              
              cat("Cumulative Mean Beta (Non-Smooth Estimation) at", sim, "iterations:","\n")
              print(round(colMeans(beta_non_smooth[1:sim, ], na.rm = TRUE), 4))
            }
            
            # Store Individual Estimates
            beta_estimates_df <- rbind(beta_estimates_df, data.frame(
              Iteration = sim,
              Tau = tau,
              Cluster_Size = cluster_size,
              Rho = rho,
              Lambda = lambda_value,
              Target_Rate = format(round(rate,2), nsmall = 2),
              Realized_Rate = format(round(realized_rates[sim],2), nsmall = 2),
              qris_Beta0 = beta_qris[sim, 1],
              qris_Beta1 = beta_qris[sim, 2],
              qris_Beta2 = beta_qris[sim, 3],
              NS_Beta0 = beta_non_smooth[sim, 1],
              NS_Beta1 = beta_non_smooth[sim, 2],
              NS_Beta2 = beta_non_smooth[sim, 3],
              CRQ_Beta0 = beta_crq[sim, 1],
              CRQ_Beta1 = beta_crq[sim, 2],
              CRQ_Beta2 = beta_crq[sim, 3],
              Xi = xi_val
            ))
          } # end of nsim simulation runs
          
          avg_qris_beta <- colMeans(beta_qris, na.rm = TRUE)
          avg_crq_beta <- colMeans(beta_crq, na.rm = TRUE)
          avg_non_smooth_beta <- colMeans(beta_non_smooth, na.rm = TRUE)
          
          conv_rate_qris <- convergence_count_qris / nsim
          conv_rate_crq <- convergence_count_crq / nsim
          conv_rate_non_smooth <- convergence_count_non_smooth / nsim
          
          results_df <- rbind(results_df, data.frame(
            Tau = tau,
            Cluster_Size = cluster_size,
            Rho = rho,
            Lambda = lambda_value,
            Target_Rate = format(round(rate, 1), nsmall = 1),
            Realized_Rate = format(round(mean(realized_rates, na.rm = TRUE), 2), nsmall = 2),
            qris_Beta0 = format(round(avg_qris_beta[1], 4), nsmall = 4, scientific = F),
            qris_Beta1 = format(round(avg_qris_beta[2], 4), nsmall = 4, scientific = F),
            qris_Beta2 = format(round(avg_qris_beta[3], 4), nsmall = 4, scientific = F),
            NS_Beta0 = format(round(avg_non_smooth_beta[1], 4), nsmall = 4, scientific = F),
            NS_Beta1 = format(round(avg_non_smooth_beta[2], 4), nsmall = 4, scientific = F),
            NS_Beta2 = format(round(avg_non_smooth_beta[3], 4), nsmall = 4, scientific = F),
            CRQ_Beta0 = format(round(avg_crq_beta[1], 4), nsmall = 4, scientific = F),
            CRQ_Beta1 = format(round(avg_crq_beta[2], 4), nsmall = 4, scientific = F),
            CRQ_Beta2 = format(round(avg_crq_beta[3], 4), nsmall = 4, scientific = F),
            Converge_qris = convergence_count_qris,
            Converge_qris_rate = format(round(conv_rate_qris, 2), nsmall = 2),
            Converge_NS = convergence_count_non_smooth,
            Converge_NS_rate = format(round(conv_rate_non_smooth, 2), nsmall = 2),
            Converge_CRQ = convergence_count_crq,
            Converge_CRQ_rate = format(round(conv_rate_crq, 2), nsmall = 2)
          ))
        } # end of targeted censoring rates
      } # end of cluster size
    } # end of rho values (intra-cluster correlations)
  } # end of taus (quantile levels)
  
  beta_estimates_df$Config <- paste(beta_estimates_df$Tau, beta_estimates_df$Rho, beta_estimates_df$Cluster_Size, round(beta_estimates_df$Lambda, 4), sep = "_")
  
  return(list(results_summary = results_df, estimates_df = beta_estimates_df))
}


##########################################################################################
# (III) Parameter settings and simulation
##########################################################################################

# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # MacbookPro
# load("data_list_wo_xi.RData")  # Load pre-generated data set
load("data_list.RData")  # Load pre-generated data set (tau=seq(0.1, 0.9, by = 0.2))
load("data_list2.RData")  # Load pre-generated data set (tau=seq(0.1, 0.9, by = 0.1))
load("data_list_case2.RData")  # Load pre-generated data set (Z and G changed)
load("data_list_case3.RData")  # Load pre-generated data set (Z and G all binary)

# Parameters
set.seed(123)
M <- 1e6
beta <- c(2, 1, 0.5)
# cluster_sizes <- c(2, 3)
cluster_sizes <- c(2)
num_clusters <- 400
# rho_values <- c(0.2, 0.5, 0.8)
rho_values <- c(0.2)
taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
# taus <- c(0.5)
# target_censoring_rates <- c(0.0, 0.1, 0.3, 0.5)
target_censoring_rates <- c(0.0, 0.1, 0.3)

nsim <- 500
msg_freq <- 100

# Run simulation 
start <- Sys.time()
results_list <- compare_methods(data_list, nsim, beta, taus, rho_values, cluster_sizes, target_censoring_rates)
end <- Sys.time()
end - start # Elapsed time: about 11.84 hours ()

warnings() # Check warning messages

##########################################################################################
# (IV) Results: tables and box plots
##########################################################################################
# Since compare_with_qris now returns a list with two elements:
# - results_summary: the summary results data frame
# - estimates_df: the individual iteration estimates for plotting
str(results_list$results_summary)
compare_summary <- results_list$results_summary

compare_summary <- compare_summary %>%
  mutate(
    Target_Rate = as.numeric(Target_Rate),
    Realized_Rate = as.numeric(Realized_Rate),
    qris_Beta0 = as.numeric(qris_Beta0),
    qris_Beta1 = as.numeric(qris_Beta1),
    qris_Beta2 = as.numeric(qris_Beta2),
    NS_Beta0 = as.numeric(NS_Beta0),
    NS_Beta1 = as.numeric(NS_Beta1),
    NS_Beta2 = as.numeric(NS_Beta2),
    CRQ_Beta0 = as.numeric(CRQ_Beta0),
    CRQ_Beta1 = as.numeric(CRQ_Beta1),
    CRQ_Beta2 = as.numeric(CRQ_Beta2),
    Converge_qris_rate = as.numeric(Converge_qris_rate),
    Converge_NS_rate = as.numeric(Converge_NS_rate),
    Converge_CRQ = Converge_CRQ/2,
    Converge_CRQ_rate = as.numeric(Converge_CRQ_rate)/2)

names(compare_summary)
str(compare_summary)

# (IV-2) Summary of warning messages
if (length(warning_log) > 0) {
  warning_log2 <- do.call(rbind, warning_log) %>% as.data.frame(stringsAsFactors = FALSE)
  
  # warning_log2 <- warning_log2 %>%
  #   mutate(
  #     across(where(is.list), ~ unlist(.)),
  #     Tau = as.numeric(Tau),
  #     Rho = as.numeric(Rho),
  #     Cluster_Size = as.numeric(Cluster_Size),
  #     Target_Rate = as.numeric(Target_Rate)
  #   )
  
  print("Warnings encountered during execution:")
  
  warning_summary <- warning_log2 %>%
    group_by(Tau, Rho, Cluster_Size, Target_Rate, Method) %>%
    summarise(Message_Count = n(), .groups = "drop")
  
  str(warning_summary)
  
} else {
  warning_summary <- tibble(
    Tau = numeric(),
    Rho = numeric(),
    Cluster_Size = numeric(),
    Target_Rate = numeric(),
    Method = character(),
    Message_Count = numeric()
  )
}

table(warning_log2$Message) 
str(warning_log2)
str(warning_summary)
head(warning_summary)
head(warning_log2)

compare_summary <- compare_summary %>%
  left_join(
    warning_summary %>%
      filter(Method == "CRQ") %>%
      rename(NonUnique_crq = Message_Count) %>%
      dplyr::select(Tau, Rho, Cluster_Size, Target_Rate, NonUnique_crq),
    by = c("Tau", "Rho", "Cluster_Size", "Target_Rate")
  ) %>%
  left_join(
    warning_summary %>%
      filter(Method == "NonSmooth") %>%
      rename(NonUnique_NS = Message_Count) %>%
      dplyr::select(Tau, Rho, Cluster_Size, Target_Rate, NonUnique_NS),
    by = c("Tau", "Rho", "Cluster_Size", "Target_Rate")
  )

compare_summary <- compare_summary %>%
  mutate(NonUnique_crq = coalesce(NonUnique_crq, 0),  
         NonUnique_crq_rate = format(round(NonUnique_crq / nsim, 2), nsmall = 2),
         NonUnique_NS = coalesce(NonUnique_NS, 0),  
         NonUnique_NS_rate = format(round(NonUnique_NS / nsim, 2), nsmall = 2)) %>%
  mutate(NonUnique_NS_rate = as.numeric(NonUnique_NS_rate),
         NonUnique_crq_rate = as.numeric(NonUnique_crq_rate))

names(compare_summary)
str(compare_summary)
head(compare_summary) 
str(warning_summary)
head(warning_summary)

write.csv(warning_summary, "warning_log_v1.2.9(wo_cens0).csv", row.names = FALSE)
write.csv(compare_summary, "results_summary_v1.2.9(wo_cens0).csv", row.names = FALSE)

# (IV-3) Boxplots of beta estimates
# Extract the individual iteration estimates for plotting
estimates_df <- results_list$estimates_df
estimates_df <- estimates_df %>%
  mutate(Target_Rate = as.numeric(Target_Rate),
         Target_Rate = round(Target_Rate, 1),
         Realized_Rate = as.numeric(Realized_Rate),
         Realized_Rate = round(Realized_Rate, 2))

head(estimates_df)
colSums(is.na(estimates_df))
unique(estimates_df$Target_Rate)


# Reshape estimates_df to long format for both QRIS and NS betas
long_estimates_df <- estimates_df %>%
  pivot_longer(cols = c(qris_Beta0, qris_Beta1, qris_Beta2,
                        NS_Beta0, NS_Beta1, NS_Beta2,
                        CRQ_Beta0, CRQ_Beta1, CRQ_Beta2), 
               names_to = c("Method", "Beta"),
               names_sep = "_",
               values_to = "Estimate") %>%
  mutate(Method = case_when(
    Method == "qris" ~ "QRIS",
    Method == "NS" ~ "NS",
    Method == "CRQ" ~ "CRQ",
    TRUE ~ Method  # This keeps any unexpected values unchanged
  ))

head(long_estimates_df)  # Check the structure

## saving diff data sets(df1=original; df2=case2; df3=case3)
# df1 <- long_estimates_df
# df2 <- long_estimates_df
# df3 <- long_estimates_df

# Compute five-number summary for each combination of Target_Rate, Method, and Tau
summary_stats <- long_estimates_df %>%
  group_by(Tau, Rho, Cluster_Size, Target_Rate, Method, Beta) %>%
  dplyr::summarise(
    Min = ifelse(all(is.na(Estimate)), NA_real_, round(min(Estimate, na.rm = TRUE),4)),
    Q1 = ifelse(all(is.na(Estimate)), NA_real_, round(quantile(Estimate, 0.25, na.rm = TRUE),4)),
    Median = ifelse(all(is.na(Estimate)), NA_real_, round(median(Estimate, na.rm = TRUE),4)),
    Q3 = ifelse(all(is.na(Estimate)), NA_real_, round(quantile(Estimate, 0.75, na.rm = TRUE),4)),
    Max = ifelse(all(is.na(Estimate)), NA_real_, round(max(Estimate, na.rm = TRUE),4)),
    .groups = "drop"
  )

str(summary_stats)

# NA counts by group
na_counts_group <- long_estimates_df %>%
  group_by(Tau, Rho, Cluster_Size, Target_Rate, Method, Beta) %>%
  summarise(na_count = sum(is.na(Estimate)), .groups = "drop") %>%
  filter(na_count > 0)

print(na_counts_group)

# Create a named vector to map Beta labels
beta_labels <- c("Beta0" = paste0("Beta0 (true = ", beta[1], ")"),
                 "Beta1" = paste0("Beta1 (true = ", beta[2], ")"),
                 "Beta2" = paste0("Beta2 (true = ", beta[3], ")"))

# Function to draw a scatter plot
gen_scp <- function(long_estimates_df, tau_val, rho_val, target_rate_val, method_sel) {
  
  # Create a data subset for the selected parameters
  scatterp_data <- long_estimates_df %>%
    filter(Tau == tau_val, Rho == rho_val, Target_Rate == target_rate_val, Method == method_sel)
  
  # Generate the scatterplot
  p <- ggplot(scatterp_data, aes(x = Iteration, y = Estimate, color = Beta)) +
    geom_point(alpha = 0.5) +
    facet_grid(Beta ~ Cluster_Size, scales = "free_y", 
               labeller = labeller(Beta = beta_labels, 
                                   Cluster_Size = c("2" = "Cluster size = 2", 
                                                    "3" = "Cluster size = 3"))) +
    theme_bw() +
    labs(title = paste0("Scatterplot of Beta Estimates over Iterations (", method_sel, 
                        ", Tau=", tau_val,
                        ", Rho=", rho_val,
                        ", Target Rate=", target_rate_val, ")"),
         x = "Iteration", y = "Beta Estimate")
  
  return(p)  # Return the plot
  
}

# Example: Generate scatterplot for QRIS method, Tau = 0.3, Target Rate = 0.1
gen_scp(long_estimates_df, tau_val = 0.3, rho_val = 0.2, target_rate_val = 0.1, method_sel = "QRIS")

gen_scp(long_estimates_df, tau_val = 0.1, rho_val = 0.2, target_rate_val = 0.5, method_sel = "QRIS")

gen_scp(long_estimates_df, tau_val = 0.3, rho_val = 0.2, target_rate_val = 0.5, method_sel = "CRQ")

gen_scp(long_estimates_df, tau_val = 0.1, rho_val = 0.2, target_rate_val = 0.1, method_sel = "NS")

gen_scp(long_estimates_df, tau_val = 0.7, rho_val = 0.2, target_rate_val = 0.5, method_sel = "NS")
gen_scp(long_estimates_df, tau_val = 0.7, rho_val = 0.5, target_rate_val = 0.5, method_sel = "NS") 
gen_scp(long_estimates_df, tau_val = 0.9, rho_val = 0.2, target_rate_val = 0.3, method_sel = "NS")
gen_scp(long_estimates_df, tau_val = 0.9, rho_val = 0.2, target_rate_val = 0.5, method_sel = "NS")

str(estimates_df)
# Function to generate box plots 
boxp_all <- function(estimates_df, beta_col, beta_value, y_label, title, cluster_size, methods_sel = c("CRQ", "NS", "QRIS")) {
  
  # Filter Data
  data_subset <- estimates_df %>%
    filter(Cluster_Size == cluster_size, Beta == beta_col, !is.na(Estimate), Method %in% methods_sel)
  
  # Define Correct Censoring Rate Order
  correct_order <- c(0, 0.1, 0.3, 0.5)
  target_censoring_rates_unique <- intersect(correct_order, unique(data_subset$Target_Rate))
  
  if (length(target_censoring_rates_unique) == 0) {
    message("No data available for Cluster Size = ", cluster_size, " and Beta = ", beta_col)
    return(NULL)
  }
  
  # Define Colors for Rho Values
  custom_colors <- c("0.2" = "skyblue", "0.5" = "yellow", "0.8" = "darkorange")
  
  # Generate Boxplots for Each Target Censoring Rate
  plots <- lapply(target_censoring_rates_unique, function(rate) {
    plot_data <- subset(data_subset, Target_Rate == rate)
    
    # # Compute Outlier Thresholds (Less Aggressive)
    # q_low <- quantile(plot_data$Estimate, 0.001, na.rm = TRUE)  # Lower 0.5%
    # q_high <- quantile(plot_data$Estimate, 0.999, na.rm = TRUE)  # Upper 0.5%
    # 
    # # Apply thresholding for visualization (but don't remove them)
    # plot_data <- plot_data %>%
    #   mutate(Extreme = ifelse(Estimate < q_low | Estimate > q_high, "Extreme", "Normal"))
    # 
    # plot_data <- plot_data %>% filter(Estimate >= q_low, Estimate <= q_high)
    
    ggplot(plot_data, aes(x = as.factor(Tau), y = Estimate, fill = as.factor(Rho))) +
      geom_boxplot(outlier.shape = 16, outlier.size = 2) +  # Keep moderate outliers
      geom_hline(yintercept = beta_value, color = "red", linewidth = 1) +
      scale_fill_manual(values = custom_colors) +  
      labs(
        x = "Tau", 
        y = y_label, 
        title = paste("Censoring Rate =", rate),
        fill = "Rho"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      ) +
      facet_grid(. ~ Method, scales = "free_y")  # Align NS & QRIS Horizontally
  })
  
  # Arrange Plots in 2x2 Grid
  grid.arrange(
    grobs = plots,
    ncol = 2,
    nrow = 2,
    top = textGrob(
      paste0(title, " (Cluster size=", cluster_size, ")"),
      gp = gpar(fontsize = 18, fontface = "bold")
    )
  )
}

# Boxplots for Beta0
boxp_all(long_estimates_df, "Beta0", beta[1], "Beta0 Estimate", "Boxplot of Beta0 estimates", cluster_size = 2)
boxp_all(long_estimates_df, "Beta0", beta[1], "Beta0 Estimate", "Boxplot of Beta0 estimates", cluster_size = 2, methods_sel = c("NS","QRIS"))

boxp_all(long_estimates_df, "Beta0", beta[1], "Beta0 Estimate", "Boxplot of Beta0 estimates", cluster_size = 3)
# Boxplots for Beta1
boxp_all(long_estimates_df, "Beta1", beta[2], "Beta1 Estimate", "Boxplot of Beta1 estimates", cluster_size = 2, methods_sel = c("NS","QRIS"))
boxp_all(long_estimates_df, "Beta1", beta[2], "Beta1 Estimate", "Boxplot of Beta1 estimates", cluster_size = 2)

boxp_all(long_estimates_df, "Beta1", beta[2], "Beta1 Estimate", "Boxplot of Beta1 estimates", cluster_size = 3)
# Boxplots for Beta2
boxp_all(long_estimates_df, "Beta2", beta[3], "Beta2 Estimate", "Boxplot of Beta2 estimates", cluster_size = 2, methods_sel = c("NS","QRIS"))
boxp_all(long_estimates_df, "Beta2", beta[3], "Beta2 Estimate", "Boxplot of Beta2 estimates", cluster_size = 3)


# Function to generate boxplots for a given method (QRIS or NS)
boxp_each <- function(estimates_df, beta_col, beta_value, y_label, title, cluster_size) {
  
  # Filter the data for the specified cluster size and beta coefficient
  data_subset <- estimates_df %>%
    filter(Cluster_Size == cluster_size, Beta == beta_col, !is.na(Estimate))
  
  # Define the correct censoring rate order
  correct_order <- c(0, 0.1, 0.3, 0.5)
  
  # Filter and ensure only the selected censoring rates are used
  data_subset <- data_subset %>%
    filter(Target_Rate %in% correct_order)
  
  # Get unique target censoring rates available in dataset, keeping order
  target_censoring_rates_unique <- intersect(correct_order, unique(data_subset$Target_Rate))
  
  # Define custom colors for Rho values
  custom_colors <- c("0.2" = "skyblue", "0.5" = "yellow", "0.8" = "darkorange")
  
  # Generate plots
  plots <- lapply(target_censoring_rates_unique, function(rate) {
    ggplot(data_subset %>% filter(Target_Rate == rate), 
           aes(x = as.factor(Tau), y = Estimate, fill = as.factor(Rho))) +
      geom_boxplot(outlier.shape = 16, outlier.size = 2) +
      geom_hline(yintercept = beta_value, color = "red", linewidth = 1) +
      scale_fill_manual(values = custom_colors) +  
      labs(
        x = "Tau", 
        y = y_label, 
        title = paste("Censoring Rate =", rate),
        fill = "Rho"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      ) +
      scale_y_continuous(limits = function(x) range(c(x, beta_value), na.rm = TRUE))
  })
  
  # Ensure exactly 4 plots (for 2x2 layout) by adding blank ones if necessary
  while (length(plots) < 4) {
    plots <- c(plots, list(ggplot() + theme_void()))  # Fix: Correctly adding empty plots
  }
  
  # Arrange the plots in a 2x2 grid with a title
  grid.arrange(
    grobs = plots,
    ncol = 2,
    nrow = 2,
    top = textGrob(
      paste0(title, " (Cluster size=", cluster_size, ")"), 
      gp = gpar(fontsize = 18, fontface = "bold")
    )
  )
}

# Draw the box plot for each function
boxp_each(subset(long_estimates_df, Method == "QRIS"), "Beta0", beta[1], "Beta0 Estimate", "Boxplot of QRIS Beta0 estimates", cluster_size = 2)
boxp_each(subset(long_estimates_df, Method == "QRIS"), "Beta0", beta[1], "Beta0 Estimate", "Boxplot of QRIS Beta0 estimates", cluster_size = 3)

boxp_each(subset(long_estimates_df, Method == "QRIS"), "Beta1", beta[2], "Beta1 Estimate", "Boxplot of QRIS Beta1 estimates", cluster_size = 2)
boxp_each(subset(long_estimates_df, Method == "QRIS"), "Beta1", beta[2], "Beta1 Estimate", "Boxplot of QRIS Beta1 estimates", cluster_size = 3)

boxp_each(subset(long_estimates_df, Method == "QRIS"), "Beta2", beta[3], "Beta2 Estimate", "Boxplot of QRIS Beta2 estimates", cluster_size = 2)
boxp_each(subset(long_estimates_df, Method == "QRIS"), "Beta2", beta[3], "Beta2 Estimate", "Boxplot of QRIS Beta2 estimates", cluster_size = 3)

boxp_each(subset(long_estimates_df, Method == "CRQ"), "Beta0", beta[1], "Beta0 Estimate", "Boxplot of CRQ Beta0 estimates", cluster_size = 2)
boxp_each(subset(long_estimates_df, Method == "CRQ"), "Beta0", beta[1], "Beta0 Estimate", "Boxplot of CRQ Beta0 estimates", cluster_size = 3)

boxp_each(subset(long_estimates_df, Method == "CRQ"), "Beta1", beta[2], "Beta1 Estimate", "Boxplot of CRQ Beta1 estimates", cluster_size = 2)
boxp_each(subset(long_estimates_df, Method == "CRQ"), "Beta1", beta[2], "Beta1 Estimate", "Boxplot of CRQ Beta1 estimates", cluster_size = 3)

boxp_each(subset(long_estimates_df, Method == "CRQ"), "Beta2", beta[3], "Beta2 Estimate", "Boxplot of CRQ Beta2 estimates", cluster_size = 2)
boxp_each(subset(long_estimates_df, Method == "CRQ"), "Beta2", beta[3], "Beta2 Estimate", "Boxplot of CRQ Beta2 estimates", cluster_size = 3)

boxp_each(subset(long_estimates_df, Method == "NS"), "Beta0", beta[1], "Beta0 Estimate", "Boxplot of NS Beta0 estimates", cluster_size = 2)
boxp_each(subset(long_estimates_df, Method == "NS"), "Beta0", beta[1], "Beta0 Estimate", "Boxplot of NS Beta0 estimates", cluster_size = 3)

boxp_each(subset(long_estimates_df, Method == "NS"), "Beta1", beta[2], "Beta1 Estimate", "Boxplot of NS Beta1 estimates", cluster_size = 2)
boxp_each(subset(long_estimates_df, Method == "NS"), "Beta1", beta[2], "Beta1 Estimate", "Boxplot of NS Beta1 estimates", cluster_size = 3)

boxp_each(subset(long_estimates_df, Method == "NS"), "Beta2", beta[3], "Beta2 Estimate", "Boxplot of NS Beta2 estimates", cluster_size = 2)
boxp_each(subset(long_estimates_df, Method == "NS"), "Beta2", beta[3], "Beta2 Estimate", "Boxplot of NS Beta2 estimates", cluster_size = 3)



# Example: before combining results
df1$Design <- "Case 1(Z:binary, G:conti.)"
df2$Design <- "Case 2(Z:conti., G:binary)"
df3$Design <- "Case 3(Z:binary, G:binary)"

# Combine all into one
long_estimates_df2 <- rbind(df1, df2, df3)

# Draw a box plot with different data sets
boxp_comb1 <- function(estimates_df, beta_col, beta_value, y_label, title, cluster_size,
                      methods_sel = c("NS", "QRIS"), save_plots = FALSE, save_dir = NULL) {
  
  tau_vals <- sort(unique(estimates_df$Tau))
  
  for (tau_val in tau_vals) {
    
    data_subset <- estimates_df %>%
      filter(Cluster_Size == cluster_size, Beta == beta_col, Tau == tau_val,
             !is.na(Estimate), Method %in% methods_sel) %>%
      mutate(
        Target_Rate = factor(Target_Rate, levels = c(0, 0.1, 0.3)),
        Rho = as.factor(Rho),
        Design = factor(
          Design,
          levels = c(
            "Case 1(Z:binary, G:conti.)",
            "Case 2(Z:conti., G:binary)",
            "Case 3(Z:binary, G:binary)"
          )
        ),
        Case_Label = case_when(
          Design == "Case 1(Z:binary, G:conti.)" ~ "Case 1",
          Design == "Case 2(Z:conti., G:binary)" ~ "Case 2",
          Design == "Case 3(Z:binary, G:binary)" ~ "Case 3"
        ),
        Rho_Method = paste0(Rho, "_", Method)
      )
    
    if (nrow(data_subset) == 0) {
      message("No data available for Tau = ", tau_val)
      next
    }
    
    custom_colors <- c(
      "0.2_NS" = "skyblue",
      "0.2_QRIS" = "steelblue",
      "0.5_NS" = "khaki1",
      "0.5_QRIS" = "goldenrod",
      "0.8_NS" = "sandybrown",
      "0.8_QRIS" = "darkorange4"
    )
    
    p <- ggplot(data_subset, aes(x = Case_Label, y = Estimate, fill = Rho_Method)) +
      geom_boxplot(outlier.shape = 16, outlier.size = 1.5) +
      geom_hline(yintercept = beta_value, color = "red", linewidth = 0.8, linetype = "dashed") +
      facet_grid(
        rows = vars(Target_Rate),
        cols = vars(Method),
        scales = "free_y",
        labeller = labeller(
          Target_Rate = c(
            "0" = "Censoring Rate = 0",
            "0.1" = "Censoring Rate = 0.1",
            "0.3" = "Censoring Rate = 0.3"
          )
        )
      ) +
      scale_fill_manual(
        name = "Rho & Method",
        values = custom_colors,
        labels = c(
          "0.2_NS" = "Rho 0.2 - NS",
          "0.2_QRIS" = "Rho 0.2 - QRIS",
          "0.5_NS" = "Rho 0.5 - NS",
          "0.5_QRIS" = "Rho 0.5 - QRIS",
          "0.8_NS" = "Rho 0.8 - NS",
          "0.8_QRIS" = "Rho 0.8 - QRIS"
        )
      ) +
      labs(
        x = "Design (Case)",
        y = y_label
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 13),
        legend.position = "right"
      )
    
    full_title <- paste0(title, " (Tau = ", tau_val, ", Cluster size = ", cluster_size, ")")
    
    if (save_plots && !is.null(save_dir)) {
      filename <- paste0(save_dir, "/plot_", beta_col, "_tau", tau_val, ".png")
      ggsave(filename, plot = p + ggtitle(full_title), width = 10, height = 6)
    } else {
      grid.arrange(
        p,
        top = textGrob(
          full_title,
          gp = gpar(fontsize = 18, fontface = "bold")
        )
      )
    }
  }
}

boxp_comb1(
  long_estimates_df2,
  "Beta0",
  beta[1],
  "Beta0 Estimate",
  "Beta0 estimates: aligned by method",
  cluster_size = 2,
  methods_sel = c("NS", "QRIS")
)


boxp_comb2 <- function(estimates_df, beta_col, beta_value, y_label, title, cluster_size,
                       methods_sel = c("NS", "QRIS"), save_plots = FALSE, save_dir = NULL) {
  
  tau_vals <- sort(unique(estimates_df$Tau))
  
  for (tau_val in tau_vals) {
    
    data_subset <- estimates_df %>%
      filter(Cluster_Size == cluster_size, Beta == beta_col, Tau == tau_val,
             !is.na(Estimate), Method %in% methods_sel) %>%
      mutate(
        Target_Rate = factor(Target_Rate, levels = c(0, 0.1, 0.3)),
        Rho = as.factor(Rho),
        Design = factor(
          Design,
          levels = c(
            "Case 1(Z:binary, G:conti.)",
            "Case 2(Z:conti., G:binary)",
            "Case 3(Z:binary, G:binary)"
          )
        ),
        Case_Label = case_when(
          Design == "Case 1(Z:binary, G:conti.)" ~ "Case 1",
          Design == "Case 2(Z:conti., G:binary)" ~ "Case 2",
          Design == "Case 3(Z:binary, G:binary)" ~ "Case 3"
        ),
        CaseMethod = factor(Case_Label, levels = c("Case 1", "Case 2", "Case 3")),
        Rho_Method = paste0(Rho, "_", Method)
      )
    
    if (nrow(data_subset) == 0) {
      message("No data available for Tau = ", tau_val)
      next
    }
    
    custom_colors <- c(
      "0.2_NS" = "skyblue",
      "0.2_QRIS" = "steelblue",
      "0.5_NS" = "khaki1",
      "0.5_QRIS" = "goldenrod",
      "0.8_NS" = "sandybrown",
      "0.8_QRIS" = "darkorange4"
    )
    
    p <- ggplot(data_subset, aes(x = CaseMethod, y = Estimate, fill = Rho_Method)) +
      geom_boxplot(
        outlier.shape = 16,
        outlier.size = 1.5,
        position = position_dodge(width = 0.8)
      ) +
      geom_hline(yintercept = beta_value, color = "red", linewidth = 0.8, linetype = "dashed") +
      facet_wrap(~ Target_Rate, nrow = 1, labeller = labeller(
        Target_Rate = c(
          "0" = "Censoring Rate = 0",
          "0.1" = "Censoring Rate = 0.1",
          "0.3" = "Censoring Rate = 0.3"
        )
      )) +
      scale_fill_manual(
        name = "Rho & Method",
        values = custom_colors,
        labels = c(
          "0.2_NS" = "Rho 0.2 - NS",
          "0.2_QRIS" = "Rho 0.2 - QRIS",
          "0.5_NS" = "Rho 0.5 - NS",
          "0.5_QRIS" = "Rho 0.5 - QRIS",
          "0.8_NS" = "Rho 0.8 - NS",
          "0.8_QRIS" = "Rho 0.8 - QRIS"
        )
      ) +
      labs(
        x = "Design (Case)",
        y = y_label
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 13),
        legend.position = "right"
      )
    
    full_title <- paste0(title, " (Tau = ", tau_val, ", Cluster size = ", cluster_size, ")")
    
    if (save_plots && !is.null(save_dir)) {
      filename <- paste0(save_dir, "/plot_", beta_col, "_tau", tau_val, ".png")
      ggsave(filename, plot = p + ggtitle(full_title), width = 12, height = 6.5)
    } else {
      grid.arrange(
        p + ggtitle(full_title),
        top = NULL
      )
    }
  }
}

boxp_comb2(
  long_estimates_df2,
  "Beta0",
  beta[1],
  "Beta0 Estimate",
  "Beta0 estimates",
  cluster_size = 2,
  methods_sel = c("NS", "QRIS")
)

boxp_comb2(
  long_estimates_df2,
  "Beta1",
  beta[2],
  "Beta1 Estimate",
  "Beta1 estimates",
  cluster_size = 2,
  methods_sel = c("NS", "QRIS")
)

boxp_comb2(
  long_estimates_df2,
  "Beta2",
  beta[3],
  "Beta2 Estimate",
  "Beta2 estimates",
  cluster_size = 2,
  methods_sel = c("NS", "QRIS")
)
