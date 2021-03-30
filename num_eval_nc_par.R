library(MASS)
library(microbenchmark)
library(foreach)
library(doParallel)
# setup for 100 replication.
nrep <- 100

# sample size
n <- 100

# will test 9 latent r and 11 zero proportion values.
latentRseq <- seq(0.05, 0.91, length.out = 9)
zratioseq <- c(0.04, 0.16, 0.28, 0.36, 0.44, 0.5, 0.56, 0.64, 0.72, 0.84, 0.96)

##### check BC
type1 <- "ternary"; type2 <- "continuous"
typesh <- "NC"
# the computation results will be saved in data.frame format
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
NC_eval <-
foreach (trueR = 1:length(latentRseq), .combine = rbind) %:%
  foreach (zrate1 = 1:length(zratioseq1), .combine = rbind){
    df <- NULL
    zratioseq2 = seq(zratioseq1[zrate1] + 0.04, 0.96, by = 0.04)
    for (zrate2 in 1:length(zratioseq2)){
      # initialize for every combination
      Kcor_org <- Kcor_ml <- Kcor_mlbd <- rep(NA, nrep)
      time_org <- time_ml <- time_mlbd <- rep(NA, nrep)
      time_all <- matrix(NA, nrow = nrep, ncol = 3)
      set.seed(123)
      for(i in 1:nrep){
        # generate bivariate normal
        z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, latentRseq[trueR], latentRseq[trueR], 1), nrow=2))
        z1shift <- quantile(z[, 1], c(zratioseq1[zrate1], zratioseq2[zrate2]))
        z1 <- z[, 1] - z1shift[1]
        z2 <- z[, 2] - z1shift[1]
        u1 <- ifelse(z1 > z1shift[2], 2, 1)
        u1[z1 <= 0] = 0          
        u2 <- z2
        # didn't apply any transformation.
        x1 <- u1
        x2 <- u2
        capture.output( # suppress the microbenchmark result.
        time_all[i, ] <- print(microbenchmark::microbenchmark(
        Kcor_org[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "original", nu = 0, tol = 1e-6)$R12,
        Kcor_ml[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, nu = 0)$R12,
        Kcor_mlbd[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "approx", nu = 0, tol = 1e-6)$R12,
        times = 10 # tried ten times
        ), unit = "us")[, 5] # to use fixed unit: "microseconds"
        # 5th column has median value
        )
      }
      # apply(time_all, 2, summary) # 3rd row gives median value.
      df <- rbind.data.frame(df, data.frame(LatentR = latentRseq[trueR], TruncRate1 = zratioseq1[zrate1], TruncRate2 = zratioseq2[zrate2], MeanAD = c(0, mean(abs(Kcor_org - Kcor_ml)), mean(abs(Kcor_org - Kcor_mlbd))), MaxAD = c(0, max(abs(Kcor_org - Kcor_ml)), max(abs(Kcor_org - Kcor_mlbd))), medianTime = apply(time_all, 2, summary)[3, ], method = c("org", "ipol", "ipol_UB")))
    }
    NC_eval <- df
  }
stopCluster(cl)
save(NC_eval, file = "NC_eval.rda")