library(MASS)
library(microbenchmark)
library(foreach)
library(doParallel)
library(Matrix)
library(chebpol)
library(pcaPP)

load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/BB_grid.rda")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/bridge.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/KendallTau.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/fromKtoR.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/estimateR.R")



# setup for 100 replication.
nrep <- 100
# sample size
n <- 100
# will test 19 latent r and 9 zero proportion values.
latentRseq <- seq(-0.9, 0.9, by = 0.1)
zratioseq <- seq(0.1, 0.9, by = 0.1)
##### check BC
type1 <- "binary"; type2 <- "binary"
# the computation results will be saved in data.frame format
cl <- makePSOCKcluster(detectCores())
registerDoParallel(cl)
BB_eval <-
foreach (trueR = 1:length(latentRseq)) %:%
  foreach (zrate = 1:length(zratioseq), .combine = rbind) %dopar% {
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
      u1 <- ifelse(z1 > 0, 1, 0)
      u2 <- ifelse(z2 > 0, 1, 0)
      # didn't apply any transformation.
      x1 <- u1
      x2 <- u2
      time_org[i] <- median(microbenchmark::microbenchmark(Kcor_org[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "original")$R12, times = 5, unit = "ms")$time)
      time_ml[i] <- median(microbenchmark::microbenchmark(Kcor_ml[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "ml")$R12, times = 5, unit = "ms")$time)
      time_mlbd[i] <- median(microbenchmark::microbenchmark(Kcor_mlbd[i] <- estimateR_mixed(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = "approx")$R12, times= 5, unit = "ms")$time)
    }
    AE <- abs(cbind(Kcor_org - latentRseq[trueR], Kcor_ml - latentRseq[trueR], Kcor_mlbd - latentRseq[trueR], Kcor_ml - Kcor_org, Kcor_mlbd - Kcor_org))
    BB_eval <- c(median(time_org), median(time_ml), median(time_mlbd), colMeans(AE), apply(AE, 2, max))
  }
stopCluster(cl)
BB_eval_3d <- array(NA, c(length(latentRseq), length(zratioseq), 13))
for (j in 1:length(latentRseq)) {
  BB_eval_3d[j, , ] = BB_eval[[j]]
}
save(BB_eval_3d, file = "BB_eval.rda")
