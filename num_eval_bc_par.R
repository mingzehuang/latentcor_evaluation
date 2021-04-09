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
zratioseq <- c(0.1, 0.9, by = 0.1)
##### check BC
type1 <- "binary"; type2 <- "continuous"
# the computation results will be saved
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
BC_eval <-
  foreach (trueR = 1:length(latentRseq)) %:%
  foreach (zrate = 1:length(zratioseq), .combine = rbind) %dopar% {
    # initialize for every combination
    time_org <- time_ml <- time_mlbd <- Kcor_org <- Kcor_ml <- Kcor_mlbd <- AE <- rep(NA, nrep)
    set.seed(123)
    for(i in 1:nrep){
      # generate bivariate normal
      z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, latentRseq[trueR], latentRseq[trueR], 1), nrow=2))
      z1shift <- quantile(z[, 1], zratioseq[zrate]) # shifting to control the truncation levels
      z1 <- z[, 1] - z1shift
      z2 <- z[, 2] - z1shift # shift the same amount of z1.
      # truncation first
      u1 <- ifelse(z1 > 0, 1, 0)
      u2 <- z2 # since this is continuous variable
      # didn't apply any transformation.
      x1 <- u1
      x2 <- u2
      time_org[i] <- median(microbenchmark::microbenchmark(Kcor_org[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "original", nu = 0, tol = 1e-6)$R12, times = 5, unit = "us")$time)
      time_ml[i] <- median(microbenchmark::microbenchmark(Kcor_ml[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "approx", nu = 0)$R12, times = 5, unit = "us")$time)
      time_mlbd[i] <- median(microbenchmark::microbenchmark(Kcor_mlbd[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "approx", nu = 0, tol = 1e-6)$R12, times= 5, unit = "us")$time)
    }
    AE <- abs(cbind(Kcor_org - latentRseq[trueR], Kcor_ml - latentRseq[trueR], Kcor_mlbd - latentRseq[trueR], Kcor_ml - Kcor_org, Kcor_mlbd - Kcor_org))
    BC_eval <- c(median(time_org), median(time_ml), median(time_mlbd), colMeans(AE), apply(AE, 2, max))
  }
stopCluster(cl)
save(BC_eval, file = "BC_eval.rda")