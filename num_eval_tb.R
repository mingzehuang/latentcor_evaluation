library(MASS)
library(microbenchmark)
# setup for 100 replication.
nrep <- 100
# sample size
n <- 100
# will test 19 latent r and 9 zero proportion values.
latentRseq <- seq(-0.9, 0.9, by = 0.1)
zratioseq <- c(0.1, 0.9, by = 0.1)
##### check TB
type1 <- "trunc"; type2 <- "binary"
# the computation results will be saved
MedianTime <- array(NA, c(length(latentRseq), length(zratioseq), 3))
MeanAE <- MaxAE <- array(NA, c(length(latentRseq), length(zratioseq), 5))
for (trueR in 1:length(latentRseq)){
  for (zrate in 1:length(zratioseq)){
    # initialize for every combination
    time_org <- time_ml <- time_mlbd <- Kcor_org <- Kcor_ml <- Kcor_mlbd <- AE <- rep(NA, nrep)
    set.seed(123)
    for(i in 1:nrep){
      # generate bivariate normal
      z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, latentRseq[trueR], latentRseq[trueR], 1), nrow=2))
      z1shift <- quantile(z[, 1], zratioseq[zrate]) # shifting to control the truncation levels
      z2shift <- quantile(z[, 2], 0.5) # fixed the zero ratio of the second variable 
      z1 <- z[, 1] - z1shift
      z2 <- z[, 2] - z2shift 
      # truncation first
      u1 <- ifelse(z1 > 0, z1, 0)
      u2 <- ifelse(z2 > 0, 1, 0)
      # didn't apply any transformation.
      x1 <- u1
      x2 <- u2
      time_Kcor_org[i] <- median(microbenchmark(Kcor_org[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "original", nu = 0, tol = 1e-6)$R12, times = 5, unit = "us")$time)
      time_Kcor_ml[i] <- median(microbenchmark(Kcor_ml[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "approx", nu = 0)$R12, times = 5, unit = "us")$time)
      time_Kcor_mlbd[i] <- median(microbenchmark(Kcor_mlbd[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "approx", nu = 0, tol = 1e-6)$R12, times= 5, unit = "us")$time)
    }
    MedianTime[trueR, zrate, ] <- c(median(time_Kcor_org), median(time_Kcor_ml), median(time_Kcor_mlbd))
    AE <- abs(rbind(Kcor_org - latentRseq[trueR], Kcor_ml - latentRseq[trueR], Kcor_mlbd - latentRseq[trueR], Kcor_ml - Kcor_org, Kcor_mlbd - Kcor_org))
    MeanAE[trueR, zrate, ] <- rowMeans(AE)
    MaxAE[trueR, zrate, ] <- max.col(AE)
  }
}
save(MedianTime, MeanAE, MaxAE, file = "TB_eval.rda")