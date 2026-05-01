# ======================== 0. Clean Environment ============================
rm(list = ls())
gc()


# ======================= 1. Environment Setup & Libraries ===================

packages.needed <- c('foreach', 'parallel', 'doParallel')
lapply(packages.needed, library, character.only = TRUE)

source('R/myFUN.R') 
source('Simulations/sim_datagen.R')


# ======================= 2. Simulation Parameters =========================

reps <- 500      
B <- 499         
alpha <- 0.05    

N1_set <- c(20, 30, 50, 60, 100)
rho_set <- c(-0.5, 0.3, 0.5, 0.7)
sigma_y_set <- c(1.5)

set.seed(20251125)


# ========================== 3. Parallel Setup =============================

no_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(no_cores, type = "SOCK")
registerDoParallel(cl)
parallel::clusterSetRNGStream(cl,)

clusterEvalQ(cl, {
  source('R/myFUN.R') 
  source('Simulations/sim_datagen.R')
})


# ================== 4. Create Directory for Results =======================

timestamp <- format(Sys.time(), "%Y%m%d-%H%M")
save_dir <- paste0("Sim_Results_Example6_Power_", timestamp)
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE) 
  cat("--- Directory created:", save_dir, "---\n")
}    


# ===================== 5. Run Parallel Simulations ========================

cat(sprintf("\n--- Starting %d simulation replications on %d cores... ---\n", reps, no_cores))

all_results_table <- NULL

# 1st Loop: Sample Size (N1)
for (N1 in N1_set) {
  cat(sprintf("\n--- Running simulations for N1 = %d ---\n", N1))
  
  current_row_results <- c(N1_M1 = N1, N2_M2 = N1)
  start_time <- proc.time()
  
  # 2nd Loop: Correlation (rho)
  for (rho in rho_set) {
    cat(sprintf("  --> Running simulations for rho = %g\n", rho))
    
    # 3rd Loop: Sigma_y
    for (sigma_y in sigma_y_set) {
      cat(sprintf("    -> Running for sigma_y = %g\n", sigma_y))
      
      result_matrix = foreach(k = 1:reps, .combine = 'rbind') %dopar% {
        
        # Balanced and Paired setup
        N2 <- N1; M1 <- N1; M2 <- N2
        
        # DGP for Example 6
        data <- data_gen5_6(N1, N2, M1, M2, mu_y = 0, sigma_y = sigma_y, rho = rho)
        X <- pnorm(data$X)
        Y <- pnorm(data$Y)
        
        # Empirical CDF
        z <- sort(unique(c(Y, X)))
        Fx <- myECDF2(X, z)$Fy
        Fy <- myECDF2(Y, z)$Fy
        diff_FyFx <- Fy - Fx
        
        # Test statistics
        N <- N1 * N2; Nbar <- min(N1, N2)
        M <- M1 * M2; Mbar <- min(M1, M2)
        
        TN_WFH <- sqrt(Nbar) * max(diff_FyFx)
        TN_DH <- sqrt(N * M / (N + M)) * max(diff_FyFx)
        
        # Recentering Funcontions
        delta_WFH <- -0.1 * sqrt(log(log(Nbar + Mbar)))
        delta_DH <- -0.1 * sqrt(log(log(N + M)))
        
        mu_WFH <- diff_FyFx * ((sqrt(Nbar) * diff_FyFx) < delta_WFH)
        mu_DH  <- diff_FyFx * ((sqrt(N) * diff_FyFx) < delta_DH)
        
        TN_bs_WFH <- numeric(B)
        TN_bs_DH  <- numeric(B)
        
        for (j in 1:B) {
          
          # WFH Bootstrap (Correlated sample)
          sample_bs <- bs_sample_cor(N1, N2, X, Y)
          Y_bs_WFH <- sample_bs$Y_bs
          X_bs_WFH <- sample_bs$X_bs
          
          Fy_bs_WFH <- myECDF2(Y_bs_WFH, z)$Fy
          Fx_bs_WFH <- myECDF2(X_bs_WFH, z)$Fy
          DN_WFH <- (Fy_bs_WFH - Fx_bs_WFH) - diff_FyFx
          
          TN_bs_WFH[j] <- sqrt(Nbar) * max(DN_WFH + mu_WFH)
          
          # DH Bootstrap (I.I.D sample)
          Y_bs_DH <- Y[sample(1:N, N, replace = TRUE)]
          X_bs_DH <- X[sample(1:M, M, replace = TRUE)]
          
          Fy_bs_DH <- myECDF2(Y_bs_DH, z)$Fy
          Fx_bs_DH <- myECDF2(X_bs_DH, z)$Fy
          DN_DH <- (Fy_bs_DH - Fx_bs_DH) - diff_FyFx
          
          TN_bs_DH[j] <- sqrt(N * M / (N + M)) * max(DN_DH + mu_DH)
        }
        
        pvalue_WFH <- (sum(TN_bs_WFH > TN_WFH) + 1) / (B + 1)
        pvalue_DH  <- (sum(TN_bs_DH > TN_DH) + 1) / (B + 1)
        
        return(c(WFH = pvalue_WFH, DH = pvalue_DH))
      }
      
      rej_rates <- apply(result_matrix, 2, function(x) mean(x < alpha))
      
      names(rej_rates) <- paste0(names(rej_rates), sprintf("(rho=%.1f,sigma=%.1f)", rho, sigma_y))
      
      current_row_results <- c(current_row_results, rej_rates)
      
    } # End of sigma_y loop
  } # End of rho loop
  
  end_time <- proc.time()
  elapsed_time <- end_time - start_time
  cat(sprintf("  -> Total time for N1=%d: %.2f seconds.\n", N1, elapsed_time[3]))
  
  all_results_table <- rbind(all_results_table, current_row_results)
  
} # End of N1 loop


# ===================== 6. Finalize and Save Results =======================

rownames(all_results_table) <- NULL

cat("\n--- Final Results Table (Preview) ---\n")
print(all_results_table[, 1:min(ncol(all_results_table), 6)])
cat("... (Full columns saved to CSV)\n")

file_name <- paste0("Table_Example6_Power_Combined_", timestamp, ".csv")
full_path <- file.path(save_dir, file_name)
write.csv(all_results_table, file = full_path, row.names = FALSE, quote = FALSE)
cat(sprintf("\n--- All simulations complete. Final results saved to: %s\n", full_path))


# ========================= 7. Stop the Cluster ============================
stopCluster(cl)
