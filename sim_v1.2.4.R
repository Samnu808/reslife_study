rm(list = ls())

##########################################################################################
# 
# 
##########################################################################################

##########################################################################################
# (I) Load necessary libraries
##########################################################################################
library(copula)
library(MASS)
library(survival)
library(quantreg)
library(ggplot2)
library(grid)
library(gridExtra)


##########################################################################################
# (II) Functions 
##########################################################################################

# (II-1) Function to generate copula-based errors for clustered survival data
gen_copula_err <- function(cluster_size, rho, marginal = "normal") {
  
  # Function to convert Kendall's tau (τ) to Pearson's correlation (ρ)
  # The transformation follows the relationship: ρ = sin(π * τ / 2)
  tau_to_rho <- function(tau) sin(pi * tau / 2)
  
  # Convert intra-cluster correlation from Kendall’s tau (τ) to Pearson’s correlation (ρ)
  rho_gaussian <- tau_to_rho(rho)
  
  # Define the copula to generate correlated errors within clusters
  if (marginal == "normal") {
    # Gaussian copula: assumes normally distributed errors with specified correlation
    cop <- normalCopula(rho_gaussian, dim = cluster_size, dispstr = "ex")
  } else if (marginal == "t") {
    # t-Copula: allows for heavier-tailed error distributions (df = 2)
    cop <- tCopula(rho_gaussian, dim = cluster_size, dispstr = "ex", df = 2)
  } else {
    stop("Unsupported marginal distribution! Choose 'normal' or 't'.")
  }
  
  # Generate a single row of correlated uniform random variables from the copula
  u <- rCopula(1, cop)
  
  # Transform uniform samples to the desired marginal distribution
  if (marginal == "normal") {
    errors <- qnorm(u)  # Convert uniform values to standard normal distribution
  } else {
    errors <- qt(u, df = 2)  # Convert uniform values to t-distribution with df = 2
  }
  
  return(as.vector(errors))  # Return the generated vector of correlated errors
}

# (II-2) Function to generate data set
gen_data <- function(cluster_size, num_clusters, tau, rho, beta, marginal, lambda) {
  data_list <- list()
  
  # Quantile adjustment (xi: shift parameter)
  xi <- if (tau == 0.5) {
    0
  } else if (marginal == "normal") {
    -qnorm(tau)
  } else {
    -qt(tau, df = 2)
  }
  
  for (i in 1:num_clusters) {
    Z <- rbinom(cluster_size, 1, 0.5)  # Individual-level covariate
    G <- rnorm(1)  # Cluster-level covariate (shared for all failure times in a subject)
    
    # Generate correlated errors for the current cluster
    errors <- gen_copula_err(cluster_size, rho, marginal)
    
    # Add the quantile adjustment
    epsilon <- errors + xi
    
    # Generate log-transformed failure times
    logT <- beta[1] + beta[2] * Z + beta[3] * G + epsilon
    
    # Generate log of censoring times
    if (lambda == 0) {
      logC <- rep(Inf, cluster_size)  # No censoring
    } else {
      logC <- log(rexp(cluster_size, rate = lambda))
    }
    
    # Observed times and event indicators
    log_obs <- pmin(logT, logC)
    delta <- as.numeric(logT <= logC)
    
    # Enforce zero censoring rate when lambda == 0 
    if (lambda == 0) {
      delta <- rep(1, cluster_size)  # Ensure all events are observed
    } else {
      delta[which.max(log_obs)] <- 1  # Ensure the largest observed time is a failure (to make G_hat not 0 later on)
    }
    
    # Combine data for the cluster
    cluster_data <- data.frame(
      cluster = rep(i, cluster_size),
      Z = Z, 
      G = rep(G, cluster_size),
      log_obs = log_obs, 
      delta = delta, 
      logC = logC, 
      logT = logT
    )
    data_list[[i]] <- cluster_data
  }
  
  # Combine data from all clusters
  result <- do.call(rbind, data_list)
  return(result)
}

# (II-3) Function to find the lambda for censoring time generation 
target_censoring_rate <- function(target_rates, taus, rho_values, num_clusters, beta, marginal, max_iter = 100, tolerance = 0.002, initial_lambda = 0.3) {
  censoring_results <- list()
  
  for (tau in taus) {
    for (rho in rho_values) {
      # Create a key for each (tau, rho) combination
      key <- paste0("Tau_", tau, "_Rho_", rho)
      censoring_results[[key]] <- data.frame(lambda = numeric(), achieved_rate = numeric())
      
      for (target_rate in target_rates) {
        # Handle target censoring rate = 0
        if (target_rate == 0) {
          lambda <- 0  # No censoring
          data <- gen_data(cluster_size = 2, num_clusters, tau, rho, beta, marginal, lambda)
          achieved_rate <- round(mean(1 - data$delta), 4)  # Achieved censoring rate (should be 0)
          
          # Save results directly without iteration
          censoring_results[[key]] <- rbind(censoring_results[[key]], 
                                            data.frame(lambda = lambda, achieved_rate = achieved_rate))
          cat("\nTarget Rate is 0: Setting Lambda to 0 | Tau:", tau, "| Rho:", rho, "| Achieved Rate:", achieved_rate, "\n")
          next  # Skip further processing for the target rate = 0
        }
        
        # Default logic for non-zero target rates
        lambda <- initial_lambda  # Initial guess
        prev_lambda <- lambda
        iteration_converged <- TRUE
        oscillation_detected <- FALSE
        
        for (iter in 1:max_iter) {
          # Generate data for a fixed cluster size(=2)
          data <- gen_data(cluster_size = 2, num_clusters, tau, rho, beta, marginal, lambda)
          achieved_rate <- round(mean(1 - data$delta), 4)
          
          # Check stopping criteria
          if (abs(achieved_rate - target_rate) < tolerance) {
            break
          }
          
          # Update lambda proportionally to the error
          adjustment <- (achieved_rate - target_rate) / target_rate
          prev_lambda <- lambda
          lambda <- lambda * (1 - adjustment)
          
          # Dynamically adjust bounds based on the target_rate
          lambda_upper_bound <- max(2, target_rate * 10)  # Empirically determined bound
          lambda <- max(min(lambda, lambda_upper_bound), 1e-4)  # Apply bounds
          
          # Check for oscillation or instability
          if (abs(lambda - prev_lambda) < 1e-6) {
            oscillation_detected <- TRUE
            iteration_converged <- FALSE
            break
          }
        }
        
        # Print only the last iteration before convergence or early stop
        cat(sprintf("\nTarget Rate is %.1f | Lambda: %.5f | Tau: %.1f | Rho: %.1f | Achieved Rate: %.1f (Iterations: %d)", 
                    target_rate, lambda, tau, rho, achieved_rate, iter), "\n")
        
        # Immediate warnings for issues
        if (oscillation_detected) {
          warning(paste("Lambda oscillation detected for Tau:", tau, ", Rho:", rho, ", Target Rate:", target_rate, ". Stopping early."))
        }
        
        if (!iteration_converged) {
          warning(paste("Iteration did not converge for Tau:", tau, ", Rho:", rho, ", Target Rate:", target_rate))
        }
        
        # Save results for this target rate
        censoring_results[[key]] <- rbind(censoring_results[[key]], 
                                          data.frame(lambda = round(lambda, 4), achieved_rate = achieved_rate))
      }
    }
  }
  
  return(censoring_results)
}


# (II-4) Non-smooth estimation function with accurate non-uniqueness tracking
non_smooth_est <- function(taus, nsim, cluster_sizes, num_clusters, rho_values, censoring_results, beta, marginal, target_censoring_rates) {
  results <- list()  # Store summary results
  all_estimates <- list()  # Store all iteration-specific results
  
  for (tau in taus) {
    for (rho in rho_values) {
      for (cluster_size in cluster_sizes) {
        # Construct the key for accessing censoring rates
        censoring_lambda <- paste0("Tau_", tau, "_Rho_", rho)
        
        if (!censoring_lambda %in% names(censoring_results)) {
          warning("Censoring rates key not found:", censoring_lambda)
          next
        }
        
        lambda_values <- censoring_results[[censoring_lambda]]$lambda
        
        if (!censoring_lambda %in% names(target_censoring_rates)) {
          stop(paste("Key not found in target_censoring_rates:", censoring_lambda))
        }
        target_rates <- target_censoring_rates[[censoring_lambda]]
        
        for (i in seq_along(lambda_values)) {
          lambda <- lambda_values[i]
          target_rate <- target_rates[i]  # Match the target rate with lambda
          
          # Initialize storage for this configuration
          estimates <- matrix(NA, nrow = nsim, ncol = length(beta))
          xi_values <- numeric(nsim)
          realized_censor_rates <- numeric(nsim)
          convergence_status <- character(nsim)  # Track convergence status
          non_uniqueness_count <- 0  # Track non-uniqueness warnings
          
          for (sim in 1:nsim) {
            data <- gen_data(cluster_size, num_clusters, tau, rho, beta, marginal, lambda)
            logT <- data$log_obs
            delta <- data$delta
            Z <- data$Z
            G <- data$G
            
            # Calculate realized censoring rate
            realized_censor_rates[sim] <- mean(1 - delta)
            
            # Calculate G_hat for survival probabilities of censoring times
            if (all(delta == 1)) {
              G_hat <- rep(1, length(logT))
            } else {
              km_fit <- survfit(Surv(logT, 1 - delta) ~ 1)
              G_hat <- stepfun(km_fit$time, c(1, km_fit$surv))(logT)
            }
            
            # Construct the design matrix U_mat and weights
            U_mat <- cbind(1, Z, G)
            weights <- delta / pmax(G_hat, 1e-6)  # Avoid G_hat values being too small
            
            # Calculate pseudo-responses
            pseudo1 <- -apply(U_mat * weights, 2, sum)
            pseudo2 <- 2 * apply(U_mat * weights * tau, 2, sum)
            
            # Construct the regression matrix and response vector
            U_reg <- rbind(U_mat, pseudo1, pseudo2)
            Y_reg <- c(logT, M, M)  # Keep the original formulation
            wt_reg <- c(weights, 1, 1)
            
            # Perform quantile regression using L1 minimization
            gamma_fit <- NULL
            withCallingHandlers({
              gamma_fit <- tryCatch({
                fit <- rq.wfit(U_reg, Y_reg, weights = wt_reg)
                convergence_status[sim] <- "Converged"
                fit$coefficients
              }, error = function(e) {
                warning("rq.wfit failed:", conditionMessage(e))
                convergence_status[sim] <- "Failed"
                return(rep(NA, ncol(U_mat)))
              })
            }, warning = function(w) {
              if (grepl("nonunique", conditionMessage(w))) {
                non_uniqueness_count <- non_uniqueness_count + 1
              }
              invokeRestart("muffleWarning")
            })
            
            # Save iteration-specific results
            if (!any(is.na(gamma_fit))) {
              estimates[sim, ] <- gamma_fit
              xi_values[sim] <- ifelse(tau == 0.5, 0, 
                                       ifelse(marginal == "normal", -qnorm(tau), -qt(tau, df = 2)))
            }
          }
          
          # Save iteration-specific estimates for visualization
          config_key <- paste(tau, rho, cluster_size, lambda, sep = "_")
          all_estimates[[config_key]] <- data.frame(
            Iteration = 1:nsim,
            Beta0 = estimates[, 1],
            Beta1 = estimates[, 2],
            Beta2 = estimates[, 3],
            Xi = xi_values,
            Censoring_Rate = realized_censor_rates
          )
          
          # Save summary results
          key <- paste("Tau:", tau, "Rho:", rho, "Cluster Size:", cluster_size, "Lambda:", lambda)
          results[[key]] <- list(
            coefficients = colMeans(estimates, na.rm = TRUE),
            xi = mean(xi_values, na.rm = TRUE),
            avg_realized_censor_rate = mean(realized_censor_rates, na.rm = TRUE),
            convergence_summary = table(convergence_status),
            non_uniqueness_count = non_uniqueness_count
          )
          
          # Print convergence and non-uniqueness summary for this iteration
          cat("\nConvergence Summary for", key, ":\n")
          print(table(convergence_status))
          cat("\nNon-uniqueness Warnings for", key, ":", non_uniqueness_count, "\n")
        }
      }
    }
  }
  
  # Combine all iteration-specific estimates into one data frame for plotting
  estimates_df <- do.call(rbind, lapply(names(all_estimates), function(key) {
    cbind(
      Config = key,
      Tau = as.numeric(strsplit(key, "_")[[1]][1]),
      Rho = as.numeric(strsplit(key, "_")[[1]][2]),
      Cluster_Size = as.numeric(strsplit(key, "_")[[1]][3]),
      Lambda = as.numeric(strsplit(key, "_")[[1]][4]),
      all_estimates[[key]]
    )
  }))
  
  return(list(results_summary = results, estimates_df = estimates_df))
}


##########################################################################################
# (III) Simulation 
##########################################################################################

# Simulation parameters
set.seed(123)
beta <- c(2, 1, 0.5)
M <- 1e6 

taus <- c(0.1,0.3,0.5,0.7,0.9)
nsim <- 500
cluster_sizes <- c(2, 3)
num_clusters <- 100
rho_values <- c(0.2,0.5,0.8)
target_censoring_rates <- c(0, 0.1, 0.3, 0.5)

# Lambda values for the target censoring rates
lambda_values <- target_censoring_rate(target_censoring_rates, taus, rho_values, num_clusters, beta, "normal")
print(lambda_values)

target_censoring_rates <- lapply(lambda_values, function(df) df$achieved_rate) # Extract target rates for non_smooth_est function
names(target_censoring_rates) <- names(lambda_values)

# Run the simulation
start_time <- Sys.time()
results_list <- non_smooth_est(taus, nsim, cluster_sizes, num_clusters, rho_values, lambda_values, beta, "normal", target_censoring_rates)
end_time <- Sys.time()
end_time - start_time # Elapsed time


##########################################################################################
# (IV) Results: summary and plots
##########################################################################################

# (IV-1) Table for beta estimates
results_summary <- results_list$results_summary # Access the `results_summary` element within `results_list`
output <- data.frame() # Initialize an empty data frame to store the results

for (key in names(results_summary)) { # Iterate over each key in `results_summary`
  # Extract parameters using a regex
  matches <- regmatches(key, regexec("Tau: ([0-9.]+) Rho: ([0-9.]+) Cluster Size: ([0-9]+) Lambda: ([0-9.eE+-]+)", key))
  
  # Check if the match is valid
  if (length(matches[[1]]) == 5) {
    tau <- as.numeric(matches[[1]][2])
    rho <- as.numeric(matches[[1]][3])
    cluster_size <- as.numeric(matches[[1]][4])
    lambda <- as.numeric(matches[[1]][5])
    
    # Retrieve coefficients and metrics
    coef <- results_summary[[key]]$coefficients
    xi <- results_summary[[key]]$xi
    avg_realized_censor_rate <- results_summary[[key]]$avg_realized_censor_rate
    
    # Append to output if coefficients exist
    if (!is.null(coef)) {
      output <- rbind(output, data.frame(
        Tau = tau,
        Cluster_Size = cluster_size,
        Rho = rho,
        Lambda = format(round(lambda, 4),nsmall=4,scientific=F),
        Censoring_Rate = format(round(avg_realized_censor_rate, 4),nsmall=4,scientific=F),
        Beta0 = format(round(coef[1], 4),nsmall=4,scientific=F),
        Beta1 = format(round(coef[2], 4),nsmall=4,scientific=F),
        Beta2 = format(round(coef[3], 4),nsmall=4,scientific=F),
        Xi = format(round(xi, 4),nsmall=4,scientific=F)
      ))
    }
  }
}

print(output) # Check the resulting data frame

# Write results to CSV
setwd("/Users/imac/Dropbox/JSLim/Research2-1/sim_results") # imac
setwd("/Users/jisunlim/Library/CloudStorage/Dropbox/JSLim/Research2-1/sim_results") # MacbookPro

write.csv(output, "results_summary_v1.2.4.csv", row.names = FALSE)

# (IV-2) Boxplots
estimates_df <- data.frame()

for (key in names(results_list$results_summary)) { # Iterate over the results_summary in results_list
  # Extract parameters from the key
  params <- unlist(strsplit(key, " "))
  tau <- as.numeric(params[2])
  rho <- as.numeric(params[4])
  cluster_size <- as.numeric(params[7])
  lambda <- round(as.numeric(gsub(",", ".", params[9])), 4)  # Convert to numeric and round
  
  # Transform the key to match the Config column format in estimates_df
  config_key <- paste(tau, rho, cluster_size, lambda, sep = "_")
  
  # Retrieve target censoring rate using the config key
  censoring_lambda_key <- paste0("Tau_", tau, "_Rho_", rho)
  if (censoring_lambda_key %in% names(target_censoring_rates)) {
    target_censor_rate <- target_censoring_rates[[censoring_lambda_key]]
  } else {
    target_censor_rate <- NA  # Set to NA if the key is not found
  }
  
  # Retrieve all individual estimates for this key
  individual_estimates <- results_list$estimates_df[results_list$estimates_df$Config == config_key, c("Beta0", "Beta1", "Beta2")]
  xi_values <- round(results_list$estimates_df[results_list$estimates_df$Config == config_key, "Xi"], 4)
  censoring_rates <- round(results_list$estimates_df[results_list$estimates_df$Config == config_key, "Censoring_Rate"], 4)
  
  # Check if the data exists
  if (nrow(individual_estimates) > 0) {
    # Add all individual estimates to the data frame
    for (i in 1:nrow(individual_estimates)) {
      estimates_df <- rbind(estimates_df, 
                            data.frame(
                              Iteration = i,
                              Tau = tau,
                              Cluster_Size = cluster_size,
                              Rho = rho,
                              Lambda = lambda,
                              Censoring_Rate = censoring_rates[i],
                              Target_Rate = target_censor_rate[which(round(lambda_values[[censoring_lambda_key]]$lambda, 4) == lambda)],
                              Beta0 = round(individual_estimates[i, "Beta0"], 4),
                              Beta1 = round(individual_estimates[i, "Beta1"], 4),
                              Beta2 = round(individual_estimates[i, "Beta2"], 4),
                              Xi = xi_values[i]
                            ))
    }
  } else {
    warning(paste("No data found for key:", key))
  }
}

head(estimates_df) # Check the resulting data frame\
str(estimates_df) # Check the structure of the updated estimates_df

# Function to generate a boxplot for a given beta parameter
gen_boxplots <- function(estimates_df, beta_col, beta_value, y_label, title, cluster_size) {
  
  # Filter the data for the specified cluster size
  data_subset <- subset(estimates_df, Cluster_Size == cluster_size)
  
  # Unique target censoring rates
  target_censoring_rates <- unique(data_subset$Target_Rate)
  
  # Define the custom colors for Rho
  custom_colors <- c("0.2" = "skyblue", "0.5" = "yellow", "0.8" = "darkorange")
  
  # Create individual plots for each target censoring rate
  plots <- lapply(target_censoring_rates, function(rate) {
    ggplot(subset(data_subset, Target_Rate == rate), 
           aes(x = as.factor(Tau), y = !!sym(beta_col), fill = as.factor(Rho))) +
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
  
  # Arrange the plots in a 2x2 grid with a corrected title
  grid.arrange(
    grobs = plots,
    ncol = 2,
    nrow = 2,
    top = textGrob(
      paste0(title, " (Cluster size=", cluster_size, ")"),  # Corrected spacing
      gp = gpar(fontsize = 18, fontface = "bold")
    )
  )
}

# Create each box plot for beta estimates
gen_boxplots(estimates_df, "Beta0", beta[1], "Beta0 Estimate", "Boxplot of Beta0 estimates", cluster_size = 2)
gen_boxplots(estimates_df, "Beta0", beta[1], "Beta0 Estimate", "Boxplot of Beta0 estimates", cluster_size = 3)

gen_boxplots(estimates_df, "Beta1", beta[2], "Beta1 Estimate", "Boxplot of Beta1 estimates", cluster_size = 2)
gen_boxplots(estimates_df, "Beta1", beta[2], "Beta1 Estimate", "Boxplot of Beta1 estimates", cluster_size = 3)

gen_boxplots(estimates_df, "Beta2", beta[3], "Beta2 Estimate", "Boxplot of Beta2 estimates", cluster_size = 2)
gen_boxplots(estimates_df, "Beta2", beta[3], "Beta2 Estimate", "Boxplot of Beta2 estimates", cluster_size = 3)


# # Function to generate a boxplot for a given beta parameter
# generate_boxplot <- function(estimates_df, beta_col, beta_value, y_label, title) {
#   ggplot(estimates_df, aes(x = as.factor(Tau), y = !!sym(beta_col), fill = as.factor(Rho))) +
#     geom_boxplot(outlier.shape = 16, outlier.size = 2) +
#     facet_wrap(~Cluster_Size, labeller = labeller(Cluster_Size = function(x) paste("Cluster size =", x))) +
#     geom_hline(yintercept = beta_value, color = "red", size = 1) +
#     labs(x = "Tau", y = y_label, title = title, fill = "Rho") +
#     theme_minimal(base_size = 16) +
#     theme(
#       strip.text = element_text(size = 18, face = "bold"),
#       legend.title = element_text(size = 16),
#       legend.text = element_text(size = 14)
#     ) +
#     scale_y_continuous(limits = function(x) {
#       range(c(x, beta_value), na.rm = TRUE)  # Adjust y-limits to include the true beta value
#     })
# }
# 
# # Generate boxplots for Beta0, Beta1, and Beta2
# plot_beta0 <- generate_boxplot(estimates_df, "Beta0", beta[1], "Beta0 Estimate", "Boxplot of Beta0 Estimates")
# plot_beta1 <- generate_boxplot(estimates_df, "Beta1", beta[2], "Beta1 Estimate", "Boxplot of Beta1 Estimates")
# plot_beta2 <- generate_boxplot(estimates_df, "Beta2", beta[3], "Beta2 Estimate", "Boxplot of Beta2 Estimates")
# 
# # Display the plots
# print(plot_beta0)
# print(plot_beta1)
# print(plot_beta2)
