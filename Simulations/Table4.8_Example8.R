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
l_set <- c(0.1, 0.5, 0.8)       

set.seed(202405) 


# ========================== 3. Parallel Setup =============================

no_cores <- detectCores() - 1

cl <- parallel::makeCluster(no_cores, type = "SOCK")
registerDoParallel(cl)

parallel::clusterSetRNGStream(cl, )

clusterEvalQ(cl, {
  source('R/myFUN.R') 
  source('Simulations/sim_datagen.R')
})


# ================== 4. Create Directory for Results =======================

timestamp <- format(Sys.time(), "%Y%m%d-%H%M")
save_dir <- paste0("Sim_Results_Example8_Size_", timestamp)
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE) 
  cat("--- Directory created:", save_dir, "---\n")
}  


# ===================== 5. Run Parallel Simulations ========================

cat(sprintf("\n--- Starting %d simulation replications on %d cores... ---\n", reps, no_cores))

all_results_table <- NULL

for(N1 in N1_set) {
  
  cat(sprintf("\n--- Running simulations for N1 = M1 = N2 = M2 = %d ---\n", N1))
  
  current_N1_results_row <- c(N1_M1 = N1, N2_M2 = N1)
  
  start_time <- proc.time()
  
  for (l in l_set) {
    cat(sprintf("  --> Running for l = %.1f\n", l))
    result_matrix = foreach(k = 1:reps, .combine = 'rbind') %dopar% {
      
      N2 <- N1; M1 <- N1; M2 <- N1
      N <- N1 * N2; M <- M1 * M2
      Nbar <- min(N1, N2); Mbar <- min(M1, M2)
      
      #DGP for Example 8
      u1<-runif(N)
      u2<-runif(N)
      X<-(u1<=l)*u1^2/l+(u1>l)*u1
      Y<-(u2<=l)*u2+(u2>l)*(l+(u2-l)^2/(1-l))
      
      # Calculate empirical CDFs
      z <- sort(unique(c(Y, X)))
      Fx <- myECDF2(X, z)$Fy 
      Fy <- myECDF2(Y, z)$Fy
      
      diff_FyFx <- Fy - Fx
      max_diff <- max(diff_FyFx)
      
      # Test statistics
      Rate_WFH <- sqrt(Nbar * Mbar / (Nbar + Mbar))
      Rate_DH <- sqrt(N * M / (N + M))
      
      TN_WFH <- Rate_WFH * max_diff
      TN_DH <- Rate_DH * max_diff
      
      # Recentering Function
      delta_WFH <- -0.1 * sqrt(log(log(Nbar + Mbar)))
      mu_WFH <- diff_FyFx * (Rate_WFH * diff_FyFx < delta_WFH)
      
      delta_DH <- -0.1 * sqrt(log(log(N + M)))
      mu_DH <- diff_FyFx * (Rate_DH * diff_FyFx < delta_DH)
      
      
      # Bootstrap simulation
      TN_bs_WFH <- numeric(B)
      TN_bs_DH <- numeric(B)
      
      for (j in 1:B) {
        # WFH bootstrap sample
        Y_bs_WFH <- bs_sample(N1, N2, Y) 
        X_bs_WFH <- bs_sample(M1, M2, X)
        DN_WFH <- myECDF2(Y_bs_WFH, z)$Fy - myECDF2(X_bs_WFH, z)$Fy - diff_FyFx
        TN_bs_WFH[j] <- Rate_WFH * max(DN_WFH+mu_WFH)
        
        # DH bootstrap sample (standard i.i.d. resampling)
        Y_bs_DH <- Y[sample(1:N, N, replace = TRUE)]
        X_bs_DH <- X[sample(1:M, M, replace = TRUE)]
        DN_DH <- myECDF2(Y_bs_DH, z)$Fy - myECDF2(X_bs_DH, z)$Fy - diff_FyFx
        TN_bs_DH[j] <- Rate_DH * max(DN_DH+mu_DH)
      }
      
      p_val_WFH <- (sum(TN_bs_WFH > TN_WFH) + 1) / (B + 1)
      p_val_DH  <- (sum(TN_bs_DH > TN_DH) + 1) / (B + 1) 
      
      return(c(WFH = p_val_WFH, DH = p_val_DH))
    }
    
    rejRate <- apply(result_matrix, 2, function(x) mean(x < alpha))
    current_N1_results_row <- c(current_N1_results_row, rejRate)
    
  }
  
  end_time <- proc.time()
  elapsed_time <- end_time - start_time
  cat(sprintf("  -> Total simulation time for N1=%d: %.2f seconds.\n", N1, elapsed_time[3]))
  
  all_results_table <- rbind(all_results_table, current_N1_results_row)
  
} 


# ===================== 6. Finalize and Save Results =======================

colnames(all_results_table) <- c("N1(M1)", "N2(M2)", 
                                 "WFH_l=0.1", "DH_l=0.1",
                                 "WFH_l=0.5", "DH_l=0.5",
                                 "WFH_l=0.8", "DH_l=0.8")

print(all_results_table)

# Save the final results table to a CSV file
file_name <- paste0("Table4-8_Size_Results_", timestamp, ".csv")
full_path <- file.path(save_dir, file_name)
write.csv(all_results_table, file = full_path, row.names = FALSE, quote = FALSE)
cat(sprintf("\n--- All simulations complete. Final results saved to: %s\n", full_path))
