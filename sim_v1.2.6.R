rm(list = ls())
##########################################################################################
# version 1.2.6
# data: Feb 13, 2025
# modifications from earlier version:
# (1) code for generating data sets were separated from here and moved to next version 
# (2) function to compare the estimation results from qris(package name:qris) is introduced 
##########################################################################################

##########################################################################################
# (I) Load necessary libraries and pre-generated data set
##########################################################################################
# packages to load
library(survival)
library(quantreg)
library(qris)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

# Working directories
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
# setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # MacbookPro

load("data_list.RData")  # Load pre-generated data set

##########################################################################################
# (II) Functions
# (1) capt_warnin: Function to capture warnings and store detailed information
# (2) compare_with_qris: to compare two methods(qris and Non-Smooth Estimation)
##########################################################################################
# (II-1) Function to capture warnings and store detailed information
warning_log <- list() # List to store warnings
capt_warnin <- function(tau, rho, cluster_size, rate, method, message) {
  warning_log <<- append(warning_log, list(
    list(Tau = tau, Rho = rho, Cluster_Size = cluster_size, Target_Rate = rate, Method = method, Message = message)
  ))
}

# (II-2) Function to compare two methods: QRIS and Non-Smooth Estimation
compare_with_qris <- function(data_list, nsim, beta) {
  results_df <- data.frame()   
  beta_estimates_df <- data.frame()  
  
  warning_log <<- data.frame(Tau = numeric(), Rho = numeric(), Cluster_Size = numeric(), 
                             Target_Rate = numeric(), Method = character(), Message = character(),
                             stringsAsFactors = FALSE)
  
  for (lambda_key in names(data_list)) {
    datasets <- data_list[[lambda_key]]  
    realized_rates <- numeric(nsim)
    
    convergence_count_qris <- 0
    convergence_count_non_smooth <- 0
    nonunique_count_qris <- 0
    nonunique_count_non_smooth <- 0  
    
    beta_qris <- matrix(NA, nrow = nsim, ncol = length(beta))
    beta_non_smooth <- matrix(NA, nrow = nsim, ncol = length(beta))
    
    tau <- as.numeric(gsub("Tau_([0-9.]+)_.*", "\\1", lambda_key))
    rho <- as.numeric(gsub(".*_Rho_([0-9.]+)_.*", "\\1", lambda_key))
    cluster_size <- as.numeric(gsub(".*_Cluster_([0-9]+)_.*", "\\1", lambda_key))
    rate <- as.numeric(gsub(".*_Rate_([0-9.]+)", "\\1", lambda_key))
    
    # Extract xi from dataset (it's the same for all rows)
    xi_val <- unique(datasets[[1]]$xi)[1]
    
    for (sim in 1:nsim) {
      data <- datasets[[sim]]
      data <- data[order(data$log_obs),]
      lambda_value <- unique(data$lambda)
      
      logT <- data$log_obs
      delta <- data$delta
      Z <- data$Z
      G <- data$G
      
      realized_rates[sim] <- mean(1 - delta)
      
      if (all(delta == 1)) {
        G_hat <- rep(1, length(logT))
      } else {
        km_fit <- survfit(Surv(exp(logT), 1 - delta) ~ 1) 
        G_hat <- stepfun(km_fit$time, c(1, km_fit$surv))(exp(logT)) 
      }
      
      U_mat <- cbind(1, Z, G)
      weights <- delta / pmax(G_hat, 1e-6)
      
      pseudo1 <- -apply(U_mat * weights, 2, sum)
      pseudo2 <- 2 * apply(U_mat * tau, 2, sum)
      
      U_reg <- rbind(U_mat, pseudo1, pseudo2)
      Y_reg <- c(logT, M, M)
      wt_reg <- c(weights, 1, 1)
      
      # ✅ REVERTED `fit_qris` LOGIC (Older Version Without Warnings)
      fit_qris <- tryCatch(
        {
          fit <- qris(Surv(exp(logT), delta) ~ Z + G, data = data, t0 = 0, Q = tau, nB = 0, method = "nonsmooth") 
          fit  # Directly return the fit object without modifying `coef()`
        },
        warning = function(w) {
          capt_warnin(tau, rho, cluster_size, rate, "QRIS", w$message)
          if (grepl("nonunique", w$message, ignore.case = TRUE)) {
            nonunique_count_qris <- nonunique_count_qris + 1  
          }
          return(NULL)
        },
        error = function(e) return(NULL)
      )
      
      # ✅ Non-Smooth Estimation (Same as Before)
      fit_non_smooth <- tryCatch(
        {
          rq.wfit(U_reg, Y_reg, weights = wt_reg)$coefficients
        },
        warning = function(w) {
          capt_warnin(tau, rho, cluster_size, rate, "NonSmooth", w$message)
          if (grepl("nonunique", w$message, ignore.case = TRUE)) {
            nonunique_count_non_smooth <- nonunique_count_non_smooth + 1  
          }
          return(rep(NA, length(beta)))
        },
        error = function(e) return(rep(NA, length(beta)))
      )
      
      # ✅ Store Results Only if Fit Exists
      if (!is.null(fit_qris)) {
        beta_qris[sim, ] <- coef(fit_qris)  # No unnecessary modifications to coef()
        convergence_count_qris <- convergence_count_qris + 1
      } else {
        beta_qris[sim, ] <- rep(NA, length(beta))
      }
      
      if (!is.null(fit_non_smooth) && !all(is.na(fit_non_smooth))) {
        beta_non_smooth[sim, ] <- fit_non_smooth
        convergence_count_non_smooth <- convergence_count_non_smooth + 1
      } else {
        beta_non_smooth[sim, ] <- rep(NA, length(beta))
      }
      
      # ✅ Store Individual Estimates
      beta_estimates_df <- rbind(beta_estimates_df, data.frame(
        Iteration = sim,
        Tau = tau,
        Cluster_Size = cluster_size,
        Rho = rho,
        Lambda = lambda_value,
        Target_Rate = rate,
        Realized_Rate = realized_rates[sim],
        qris_Beta0 = beta_qris[sim, 1],
        qris_Beta1 = beta_qris[sim, 2],
        qris_Beta2 = beta_qris[sim, 3],
        NS_Beta0 = beta_non_smooth[sim, 1],
        NS_Beta1 = beta_non_smooth[sim, 2],
        NS_Beta2 = beta_non_smooth[sim, 3],
        Xi = xi_val
      ))
    }
    
    # ✅ Compute Summary Statistics
    avg_qris_beta <- colMeans(beta_qris, na.rm = TRUE)
    avg_non_smooth_beta <- colMeans(beta_non_smooth, na.rm = TRUE)
    
    conv_rate_qris <- convergence_count_qris / nsim
    conv_rate_non_smooth <- convergence_count_non_smooth / nsim
    
    # ✅ Store Summary Results
    results_df <- rbind(results_df, data.frame(
      Tau = tau,
      Cluster_Size = cluster_size,
      Rho = rho,
      Lambda = lambda_value,
      Target_Rate = rate,
      Realized_Rate = round(mean(realized_rates, na.rm = TRUE), 2),
      qris_Beta0 = format(round(avg_qris_beta[1], 4), nsmall = 4, scientific = F),
      qris_Beta1 = format(round(avg_qris_beta[2], 4), nsmall = 4, scientific = F),
      qris_Beta2 = format(round(avg_qris_beta[3], 4), nsmall = 4, scientific = F),
      NS_Beta0 = format(round(avg_non_smooth_beta[1], 4), nsmall = 4, scientific = F),
      NS_Beta1 = format(round(avg_non_smooth_beta[2], 4), nsmall = 4, scientific = F),
      NS_Beta2 = format(round(avg_non_smooth_beta[3], 4), nsmall = 4, scientific = F),
      Converge_qris = convergence_count_qris,
      Converge_qris_rate = format(round(conv_rate_qris, 2), nsmall = 2),
      Converge_NS = convergence_count_non_smooth,
      Converge_NS_rate = format(round(conv_rate_non_smooth, 2), nsmall = 2),
      NonUnique_qris = nonunique_count_qris,
      NonUnique_qris_rate = format(round(nonunique_count_qris / nsim, 2), nsmall = 2),
      NonUnique_NS = nonunique_count_non_smooth,
      NonUnique_NS_rate = format(round(nonunique_count_non_smooth / nsim, 2), nsmall = 2)
    ))
  }
  
  beta_estimates_df$Config <- paste(beta_estimates_df$Tau, beta_estimates_df$Rho, beta_estimates_df$Cluster_Size, 
                                    round(beta_estimates_df$Lambda, 4), sep = "_")
  
  return(list(results_summary = results_df, estimates_df = beta_estimates_df))
}


##########################################################################################
# (III) Parameter settings and simulation
##########################################################################################

# Parameters
set.seed(123)
M <- 1e6
beta <- c(2, 1, 0.5)
taus <- c(0.1, 0.3, 0.5, 0.7, 0.9)
nsim <- 10 
cluster_sizes <- c(2, 3)
num_clusters <- 100
rho_values <- c(0.2, 0.5, 0.8)
target_censoring_rates <- c(0, 0.1, 0.3, 0.5)

# Run simulation 
start <- Sys.time()
results_list <- compare_with_qris(data_list, nsim, beta)
end <- Sys.time()
end - start # Elapsed time: about 9 mins
warnings() # Check warning messages

##########################################################################################
# (IV) Results: tables and box plots
##########################################################################################

# Since compare_with_qris now returns a list with two elements:
# - results_summary: the summary results data frame
# - estimates_df: the individual iteration estimates for plotting
print(results_list$results_summary) # Print the summary table for review
str(warning_log)

# (IV-1) Summary of warning messages
if (length(warning_log) > 0) {
  # Convert warning_log (list) into a structured data frame
  warning_log <- do.call(rbind, warning_log) %>% as.data.frame(stringsAsFactors = FALSE)
  
  # Ensure numeric columns remain numeric after unlisting
  warning_log <- warning_log %>%
    mutate(
      across(where(is.list), ~ unlist(.)),
      Tau = as.numeric(Tau),
      Rho = as.numeric(Rho),
      Cluster_Size = as.numeric(Cluster_Size),
      Target_Rate = as.numeric(Target_Rate)
    )
  
  print("Warnings encountered during execution:")
  print(str(warning_log))
  
  # Summarize warning counts by key variables
  warning_summary <- warning_log %>%
    group_by(Tau, Rho, Cluster_Size, Target_Rate, Method) %>%
    summarise(Message_Count = n(), .groups = "drop")  # Ensuring no NA values
  
  print(warning_summary)
  
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

table(warning_log$Message) 
str(warning_log)
str(warning_summary)
head(warning_summary)
head(warning_log)

# (IV-2) Table of results: merged with warnings
# Make sure compare_methods_summary exists by assigning the summary results
compare_summary <- results_list$results_summary

# Merge Non-Uniqueness counts into results_summary
compare_summary <- results_list$results_summary %>%
  select(-matches("NonUnique_qris|NonUnique_NS|Method")) %>%  # Remove any pre-existing non-unique columns
  left_join(
    warning_summary %>% filter(Method == "QRIS") %>% 
      rename(NonUnique_qris = Message_Count) %>% 
      select(-Method),  # Drop the Method column
    by = c("Tau", "Rho", "Cluster_Size", "Target_Rate")
  ) %>%
  left_join(
    warning_summary %>% filter(Method == "NonSmooth") %>% 
      rename(NonUnique_NS = Message_Count) %>% 
      select(-Method),  # Drop the Method column
    by = c("Tau", "Rho", "Cluster_Size", "Target_Rate")
  ) %>%
  mutate(
    NonUnique_qris = coalesce(NonUnique_qris, 0),  
    NonUnique_qris_rate = format(round(NonUnique_qris / nsim, 2), nsmall = 2),
    NonUnique_NS = coalesce(NonUnique_NS, 0),  
    NonUnique_NS_rate = format(round(NonUnique_NS / nsim, 2), nsmall = 2)
  ) %>%
  mutate(
    across(c(qris_Beta0, qris_Beta1, qris_Beta2, 
             NS_Beta0, NS_Beta1, NS_Beta2, 
             Converge_qris_rate, Converge_NS_rate, 
             NonUnique_qris_rate, NonUnique_NS_rate), 
           ~ as.numeric(.))
  ) %>%
  select(
    Tau, Cluster_Size, Rho, Lambda, Target_Rate, Realized_Rate, 
    qris_Beta0, qris_Beta1, qris_Beta2, 
    NS_Beta0, NS_Beta1, NS_Beta2, 
    Converge_qris, Converge_qris_rate, 
    Converge_NS, Converge_NS_rate, 
    NonUnique_qris, NonUnique_qris_rate, 
    NonUnique_NS, NonUnique_NS_rate  # Ensure correct order
  )

str(compare_summary)
head(compare_summary)

# # summary for NS only
# result_NS <- compare_summary %>%
#   select(
#     Tau, Cluster_Size, Rho, Lambda, Target_Rate, Realized_Rate, 
#     NS_Beta0, NS_Beta1, NS_Beta2,
#     Converge_NS, Converge_NS_rate,
#     NonUnique_NS, NonUnique_NS_rate
#   )

# Save to a CSV file
write.csv(warning_summary, "warning_log_both(v1.2.6).csv", row.names = FALSE)
write.csv(compare_summary, "results_summary_both(v1.2.6).csv", row.names = FALSE)
# write.csv(result_NS, "results_summary_NS(v1.2.6).csv", row.names = FALSE)

# (IV-3) Boxplots of beta estimates
# Extract the individual iteration estimates for plotting
estimates_df <- results_list$estimates_df

str(estimates_df)
head(estimates_df)
colSums(is.na(estimates_df))

# Reshape estimates_df to long format for both QRIS and NS betas
long_estimates_df <- estimates_df %>%
  pivot_longer(cols = c(qris_Beta0, qris_Beta1, qris_Beta2, NS_Beta0, NS_Beta1, NS_Beta2), 
               names_to = c("Method", "Beta"),
               names_sep = "_",
               values_to = "Estimate") %>%
  mutate(Method = ifelse(Method == "qris", "QRIS", "NS"))
head(long_estimates_df)  # Check the structure

# Function to generate a boxplot for a given beta parameter
# Function to generate a boxplot for a given beta parameter
gen_boxplots_both <- function(estimates_df, beta_col, beta_value, y_label, title, cluster_size) {
  # Filter the data for the specified cluster size and beta coefficient
  data_subset <- subset(estimates_df, Cluster_Size == cluster_size & Beta == beta_col)
  
  # Remove NA values to avoid warnings
  data_subset <- data_subset %>% filter(!is.na(Estimate))
  
  # Unique target censoring rates
  target_censoring_rates_unique <- unique(data_subset$Target_Rate)
  
  # Define custom colors for QRIS and NS methods
  custom_colors <- c("QRIS" = "skyblue", "NS" = "darkorange")
  
  # Create individual plots for each target censoring rate
  plots <- lapply(target_censoring_rates_unique, function(rate) {
    ggplot(subset(data_subset, Target_Rate == rate), 
           aes(x = as.factor(Tau), y = Estimate, fill = Method)) +
      geom_boxplot(outlier.shape = 16, outlier.size = 2, position = position_dodge(0.8)) +
      geom_hline(yintercept = beta_value, color = "red", linewidth = 1) +
      scale_fill_manual(values = custom_colors) +  # Apply custom colors
      labs(
        x = "Tau", 
        y = y_label, 
        title = paste("Censoring Rate =", rate),
        fill = "Method"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      ) +
      scale_y_continuous(limits = function(x) range(c(x, beta_value), na.rm = TRUE))
  })
  
  # If there are fewer than 4 plots, add empty plots to fill a 2x2 grid
  while (length(plots) < 4) {
    plots <- c(plots, list(ggplot() + theme_void()))
  }
  
  # Arrange the plots in a 2x2 grid with a common title
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
gen_boxplots_both(long_estimates_df, "Beta0", beta[1], "Beta0 Estimate", "Boxplot of Beta0 estimates", cluster_size = 2)
gen_boxplots_both(long_estimates_df, "Beta0", beta[1], "Beta0 Estimate", "Boxplot of Beta0 estimates", cluster_size = 3)
# Boxplots for Beta1
gen_boxplots_both(long_estimates_df, "Beta1", beta[2], "Beta1 Estimate", "Boxplot of Beta1 estimates", cluster_size = 2)
gen_boxplots_both(long_estimates_df, "Beta1", beta[2], "Beta1 Estimate", "Boxplot of Beta1 estimates", cluster_size = 3)
# Boxplots for Beta2
gen_boxplots_both(long_estimates_df, "Beta2", beta[3], "Beta2 Estimate", "Boxplot of Beta2 estimates", cluster_size = 2)
gen_boxplots_both(long_estimates_df, "Beta2", beta[3], "Beta2 Estimate", "Boxplot of Beta2 estimates", cluster_size = 3)


# Function to draw boxplot for each method
gen_boxplots <- function(estimates_df, beta_col, beta_value, y_label, title, cluster_size) {
  
  # Filter the data for the specified cluster size and beta coefficient
  data_subset <- subset(estimates_df, Cluster_Size == cluster_size & Beta == beta_col)
  
  # Remove NA values
  data_subset <- data_subset %>% filter(!is.na(Estimate))
  
  # Unique target censoring rates
  target_censoring_rates <- unique(data_subset$Target_Rate)
  
  # Define the custom colors for Rho values
  custom_colors <- c("0.2" = "skyblue", "0.5" = "yellow", "0.8" = "darkorange")
  
  # Create individual plots for each target censoring rate
  plots <- lapply(target_censoring_rates, function(rate) {
    ggplot(subset(data_subset, Target_Rate == rate), 
           aes(x = as.factor(Tau), y = Estimate, fill = as.factor(Rho))) +
      geom_boxplot(outlier.shape = 16, outlier.size = 2) +
      geom_hline(yintercept = beta_value, color = "red", linewidth = 1) +
      scale_fill_manual(values = custom_colors) +  # Apply custom colors
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
  
  # Ensure there are exactly 4 plots (for 2x2 layout)
  while (length(plots) < 4) {
    plots <- c(plots, ggplot() + theme_void())
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


gen_boxplots(subset(long_estimates_df, Method == "QRIS"), "Beta0", beta[1], "Beta0 Estimate", "Boxplot of QRIS Beta0 estimates", cluster_size = 2)
gen_boxplots(subset(long_estimates_df, Method == "QRIS"), "Beta0", beta[1], "Beta0 Estimate", "Boxplot of QRIS Beta0 estimates", cluster_size = 3)

gen_boxplots(subset(long_estimates_df, Method == "QRIS"), "Beta1", beta[2], "Beta1 Estimate", "Boxplot of QRIS Beta1 estimates", cluster_size = 2)
gen_boxplots(subset(long_estimates_df, Method == "QRIS"), "Beta1", beta[2], "Beta1 Estimate", "Boxplot of QRIS Beta1 estimates", cluster_size = 3)

gen_boxplots(subset(long_estimates_df, Method == "QRIS"), "Beta2", beta[3], "Beta2 Estimate", "Boxplot of QRIS Beta2 estimates", cluster_size = 2)
gen_boxplots(subset(long_estimates_df, Method == "QRIS"), "Beta2", beta[3], "Beta2 Estimate", "Boxplot of QRIS Beta2 estimates", cluster_size = 3)

gen_boxplots(subset(long_estimates_df, Method == "NS"), "Beta0", beta[1], "Beta0 Estimate", "Boxplot of NS Beta0 estimates", cluster_size = 2)
gen_boxplots(subset(long_estimates_df, Method == "NS"), "Beta0", beta[1], "Beta0 Estimate", "Boxplot of NS Beta0 estimates", cluster_size = 3)

gen_boxplots(subset(long_estimates_df, Method == "NS"), "Beta1", beta[2], "Beta1 Estimate", "Boxplot of NS Beta1 estimates", cluster_size = 2)
gen_boxplots(subset(long_estimates_df, Method == "NS"), "Beta1", beta[2], "Beta1 Estimate", "Boxplot of NS Beta1 estimates", cluster_size = 3)

gen_boxplots(subset(long_estimates_df, Method == "NS"), "Beta2", beta[3], "Beta2 Estimate", "Boxplot of NS Beta2 estimates", cluster_size = 2)
gen_boxplots(subset(long_estimates_df, Method == "NS"), "Beta2", beta[3], "Beta2 Estimate", "Boxplot of NS Beta2 estimates", cluster_size = 3)

