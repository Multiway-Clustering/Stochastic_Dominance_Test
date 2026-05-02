# ======================== 0. Clean Environment ============================
rm(list = ls())
gc()


# ======================= 1. Environment Setup & Libraries ===================

# Load necessary packages
packages.needed <- c('foreach', 'parallel', 'doParallel')
lapply(packages.needed, library, character.only = TRUE)

# Load core functions (please adjust paths if necessary)
source('R/myFUN.R') 
source('Simulations/sim_datagen.R')


# ======================= 2. Simulation Parameters =========================

reps <- 500      # Number of simulation replications
B <- 499         # Number of bootstrap samples
alpha <- 0.05    # Significance level

mu_y_set <- c(0.2, 0, -0.2)
sigma_y_fixed <- 1 

# Block 1: Balanced cases
scenarios_balanced <- data.frame(N1 = c(20, 30, 50, 60, 100), 
                                 N2 = c(20, 30, 50, 60, 100))
# Block 2: Unbalanced cases
scenarios_unbalanced <- data.frame(N1 = c(20, 50, 60, 100), 
                                   N2 = 30)
all_scenarios <- rbind(scenarios_balanced, scenarios_unbalanced)

set.seed(20240503)


# ========================== 3. Parallel Setup =============================


no_cores <- detectCores() - 1
cl <- parallel::makeCluster(no_cores, type = "SOCK")
registerDoParallel(cl)

parallel::clusterSetRNGStream(cl,)

clusterEvalQ(cl, {
  source('R/myFUN.R') 
  source('Simulations/sim_datagen.R')
})


# ================== 4. Create Directory for Results =======================

timestamp <- format(Sys.time(), "%Y%m%d-%H%M")
save_dir <- paste0("Sim_Results_Example1_SizePower_", timestamp)
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE) 
  cat("--- Directory created:", save_dir, "---\n")
}  


# ===================== 5. Run Parallel Simulations ========================

cat(sprintf("\n--- Starting %d simulation replications on %d cores... ---\n", reps, no_cores))
cat(sprintf("--- Total scenarios to run: %d ---\n", nrow(all_scenarios)))

# Main data frame to store all final results
all_results_table <- NULL

# Loop through each scenario (each row of the table)
for (i in 1:nrow(all_scenarios)) {
  
  N1 <- all_scenarios$N1[i]
  N2 <- all_scenarios$N2[i]
  
  cat(sprintf("\n--- Running scenario: N1=%d, N2=%d ---\n", N1, N2))
  current_scenario_results_row <- c(N1_M1 = N1, N2_M2 = N2)
  
  start_time <- proc.time()
  for (mu_y in mu_y_set) {
    cat(sprintf("  --> Running for mu_y = %.1f\n", mu_y))
    
    result_matrix = foreach(k = 1:reps, .combine = 'rbind') %dopar% {
      
      M1 <- N1; M2 <- N2
      
      # DGP for Example 1
      data <- data_gen1_2(N1, N2, M1, M2, mu_y = mu_y, sigma_y = sigma_y_fixed)
      Y<-data$Y;X<-data$X
      N <- length(Y); M <- length(X)
      Nbar <- min(N1, N2); Mbar <- min(M1, M2)
      
      # Calculate empirical CDFs
      z <- sort(unique(c(Y, X)))
      Fx <- myECDF(X, z)$Fy
      Fy <- myECDF(Y, z)$Fy
      diff_FyFx <- Fy - Fx
      
      # Test statistics
      Rate_WFH <- sqrt(Nbar * Mbar / (Nbar + Mbar))
      Rate_DH <- sqrt(N * M / (N + M))
      TN_WFH <- Rate_WFH * max(diff_FyFx)
      TN_DH <- Rate_DH * max(diff_FyFx)
      
      # Recentering Function
      delta_WFH <- -0.1 * sqrt(log(log(Nbar + Mbar)))
      delta_DH <- -0.1 * sqrt(log(log(N + M)))
      mu_WFH <- diff_FyFx * (Rate_WFH * diff_FyFx < delta_WFH)
      mu_DH <- diff_FyFx * (Rate_DH * diff_FyFx < delta_DH)
      
      # Bootstrap simulation
      TN_bs_WFH <- numeric(B)
      TN_bs_DH <- numeric(B)
      
      for (j in 1:B) {
        # WFH bootstrap uses bs_sample2
        Y_bs_WFH <- bs_sample(N1, N2, Y)
        X_bs_WFH <- bs_sample(M1, M2, X)
        DN_WFH <- myECDF2(Y_bs_WFH, z)$Fy - myECDF2(X_bs_WFH, z)$Fy - diff_FyFx
        TN_bs_WFH[j] <- Rate_WFH * max(DN_WFH+mu_WFH)
        
        # DH bootstrap uses standard i.i.d. resampling
        Y_bs_DH <- Y[sample(1:N, N, replace = TRUE)]
        X_bs_DH <- X[sample(1:M, M, replace = TRUE)]
        DN_DH <- myECDF2(Y_bs_DH, z)$Fy - myECDF2(X_bs_DH, z)$Fy - diff_FyFx
        TN_bs_DH[j] <- Rate_DH * max(DN_DH + mu_DH)
      }
      
      # Calculate p-values for the two methods in the table
      p_val_WFH <- (sum(TN_bs_WFH > TN_WFH) + 1) / (B + 1)
      p_val_DH  <- (sum(TN_bs_DH > TN_DH) + 1) / (B + 1)
      
      return(c(WFH = p_val_WFH, DH = p_val_DH))
    }
    
    # Calculate and store rejection rates for the current (N1, N2, mu_y) setting
    rejRate <- apply(result_matrix, 2, function(x) mean(x < alpha))
    current_scenario_results_row <- c(current_scenario_results_row, rejRate)
    
  } 
  
  end_time <- proc.time()
  elapsed_time <- end_time - start_time
  cat(sprintf("  -> Total time for scenario (N1=%d, N2=%d): %.2f seconds.\n", N1, N2, elapsed_time[3]))
  
  all_results_table <- rbind(all_results_table, current_scenario_results_row)
  
} 


# ===================== 6. Finalize and Save Results =======================

# Set column names for the final table to match the paper's format
colnames(all_results_table) <- c("N1(M1)", "N2(M2)", 
                                 "WFH_mu=0.2", "DH_mu=0.2",
                                 "WFH_mu=0", "DH_mu=0",
                                 "WFH_mu=-0.2", "DH_mu=-0.2")

# Display final results in the console
cat("\n--- Final Results Table ---\n")
print(all_results_table)

# Save the final results table to a CSV file
file_name <- paste0("Table4-1_SizePower_Results_", timestamp, ".csv")
full_path <- file.path(save_dir, file_name)
write.csv(all_results_table, file = full_path, row.names = FALSE, quote = FALSE)
cat(sprintf("\n--- All simulations complete. Final results saved to: %s\n", full_path))


# ========================= 7. Stop the Cluster ============================
stopCluster(cl)

# Optional: Make a sound to notify completion
# if ("beepr" %in% installed.packages()) { library(beepr); beep(3) }
