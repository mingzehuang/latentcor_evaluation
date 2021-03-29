library(MASS)
library(microbenchmark)
# setup for 100 replication.
nrep <- 100

# sample size
n <- 100

# will test 9 latent r values.
latentRseq <- seq(0.05, 0.91, length.out = 9)
zratioseq1 <- seq(0.04, 0.92, by = 0.04)

##### check two cases of NC, NB
    type1 <- "ternary"; type2 <- "continuous"
  
  # the computation results will be saved in data.frame format
  df_comptime <- df_accuracy <- NULL
  
  
  for (trueR in latentRseq){
    
    for (zrate1 in zratioseq1){
      zratioseq2 = seq(zrate1 + 0.04, 0.96, by = 0.04)
      for (zrate2 in zratioseq2){
      # initialize for every combination
      Kcor_org <- Kcor_ml <- Kcor_mlbd <- rep(NA, nrep)
      time_org <- time_ml <- time_mlbd <- rep(NA, nrep)
      time_all <- matrix(NA, nrow = nrep, ncol = 3)
      
      ptm <- proc.time()
      set.seed(123)
      
      for(i in 1:nrep){
        # generate bivariate normal
        z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, trueR, trueR, 1), nrow=2))
          z1shift <- quantile(z[, 1], c(zrate1, zrate2))
          z1 <- z[, 1] - z1shift[1]
          z2 <- z[, 2] - z1shift[1]
          
          u1 <- ifelse(z1 > z1shift[2], 2, 1)
          u1[z1 <= 0] = 0          
          u2 <- z2
        
        # didn't apply any transformation.
        x1 <- u1
        x2 <- u2
        
        
        capture.output( # suppress the microbenchmark result.
          time_all[i, ] <- print(microbenchmark(
            Kcor_org[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "original", nu = 0, tol = 1e-6)$R12,
            Kcor_ml[i] <- estimateR_mixed_mlonly(X1 = x1, X2 = x2, type1 = type1, type2 = type2, nu = 0)$R12,
            Kcor_mlbd[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "approx", nu = 0, tol = 1e-6)$R12,
            times = 10 # tried ten times
          ), unit = "us")[, 5] # to use fixed unit: "microseconds"
          # 5th column has median value
        )
      }
      
      # apply(time_all, 2, summary) # 3rd row gives median value.
      df_comptime <- rbind.data.frame(df_comptime, data.frame(LatentR = trueR, TruncRate = zrate, medianTime = apply(time_all, 2, summary)[3, ], method = c("org", "ipol", "ipol_UB")))
      
      # save two kinds of errors: maximum absolute error and mean absolute error.
      df_accuracy <- rbind.data.frame(df_accuracy, data.frame(LatentR = trueR, TruncRate = zrate, MeanAD = c(mean(abs(Kcor_org - Kcor_ml)), mean(abs(Kcor_org - Kcor_mlbd))), MaxAD = c(max(abs(Kcor_org - Kcor_ml)), max(abs(Kcor_org - Kcor_mlbd))), method = c("ipol", "ipol_UB")))
      
      cat(typesh, "case: trueR = ", trueR, "\t zrate =", zrate, "\t took ", (proc.time() - ptm)[3], " seconds.\n")
      }
    }
  }
  
  save(df_comptime, df_accuracy, file = paste0("Data/TwoSim_", typesh, "_rep10.Rda"))