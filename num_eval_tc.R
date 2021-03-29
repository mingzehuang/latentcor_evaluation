library(MASS)
library(microbenchmark)
# setup for 100 replication.
nrep <- 100

# sample size
n <- 100

# will test 9 latent r and 11 zero proportion values.
latentRseq <- seq(0.05, 0.91, length.out = 9)
zratioseq <- c(0.04, 0.16, 0.28, 0.36, 0.44, 0.5, 0.56, 0.64, 0.72, 0.84, 0.96)
##### check six cases of TC, TT, BC, BB, TB, NC
type1 <- "trunc"; type2 <- "continuous"
typesh <- "TC"
# the computation results will be saved in data.frame format
df_comptime <- df_accuracy <- NULL
for (trueR in latentRseq){
  for (zrate in zratioseq){
    # initialize for every combination
    Kcor_org <- Kcor_ml <- Kcor_mlbd <- rep(NA, nrep)
    time_org <- time_ml <- time_mlbd <- rep(NA, nrep)
    time_all <- matrix(NA, nrow = nrep, ncol = 3)
    ptm <- proc.time()
    set.seed(123)
    for(i in 1:nrep){
      # generate bivariate normal
      z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, trueR, trueR, 1), nrow=2))
      # shifting to control the truncation levels
      z1shift <- quantile(z[, 1], zrate)
      z1 <- z[, 1] - z1shift
      z2 <- z[, 2] - z1shift # shift the same amount of z1.
      # truncate the first variable.
      u1 <- ifelse(z1 > 0, z1, 0)
      u2 <- z2 # since this is continuous variable
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
save(df_comptime, df_accuracy, file = paste0("TwoSim_", typesh, "_rep10.rda"))