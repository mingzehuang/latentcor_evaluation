library(MASS)
library(microbenchmark)
library(foreach)
library(doParallel)
# setup for 100 replication.
nrep <- 100
# sample size
n <- 100
# will test 19 latent r and 9 zero proportion values.
latentRseq <- seq(-0.9, 0.9, by = 0.1)
zratioseq <- seq(0.2, 0.8, by = 0.1)
##### check NN
type1 <- "ternary"; type2 <- "ternary"
# the computation results will be saved in data.frame format
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
NN_eval <-
  foreach (trueR = 1:length(latentRseq)) %:%
    foreach (zrate = 1:length(zratioseq), .combine = rbind) %dopar% {
      zrate2 <- zratioseq[zrate] + (1 - zratioseq[zrate]) / 2
      zrate21 <- zratioseq[zrate] / 2
      zrate22 <- zrate2 / 2
      # initialize for every combination
      time_org <- time_ml <- time_mlbd <- Kcor_org <- Kcor_ml <- Kcor_mlbd <- AE <- rep(NA, nrep)          set.seed(123)
      for(i in 1:nrep){
        # generate bivariate normal
        z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, trueR, trueR, 1), nrow=2))
        z1shift <- quantile(z[, 1], c(zrate11, zrate12))
        z2shift <- quantile(z[, 2], c(zrate21, zrate22))
        z1 <- z[, 1] - z1shift[1]
        z2 <- z[, 2] - z2shift[1]
        u1 <- ifelse(z1 > z1shift[2], 2, 1)
        u1[z1 <= 0] = 0
        u2 <- ifelse(z2 > z2shift[2], 2, 1)
        u2[z2 <= 0] = 0
      # didn't apply any transformation.
      x1 <- u1
      x2 <- u2
      time_org[i] <- median(microbenchmark::microbenchmark(Kcor_org[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "original", nu = 0)$R12, times = 5, unit = "ms")$time)
      time_ml[i] <- median(microbenchmark::microbenchmark(Kcor_ml[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "ml", nu = 0)$R12, times = 5, unit = "ms")$time)
      time_mlbd[i] <- median(microbenchmark::microbenchmark(Kcor_mlbd[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "approx", nu = 0)$R12, times= 5, unit = "ms")$time)
    }
    AE <- abs(cbind(Kcor_org - latentRseq[trueR], Kcor_ml - latentRseq[trueR], Kcor_mlbd - latentRseq[trueR], Kcor_ml - Kcor_org, Kcor_mlbd - Kcor_org))
    NN_eval <- c(median(time_org), median(time_ml), median(time_mlbd), colMeans(AE), apply(AE, 2, max))
    }
stopCluster(cl)
NN_eval_3d <- array(NA, c(length(latentRseq), length(zratioseq), 13))
for (j in 1:length(latentRseq)) {
  NN_eval_3d[j, , ] = NN_eval[[j]]
}
save(NN_eval_3d, file = "NN_eval.rda")