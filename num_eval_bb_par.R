library(MASS)
library(microbenchmark)
library(foreach)
library(Matrix)
library(chebpol)
library(pcaPP)
library(doFuture)

load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/all_grid.rda")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/internal.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/latentcor.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/gen_data.R")

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
registerDoFuture()
plan(multicore, workers = 80)

BB_eval <-
foreach (trueR = 1:length(latentRseq)) %:%
  foreach (zrate = 1:length(zratioseq), .combine = rbind) %dopar% {
    # initialize for every combination
    time_org <- time_ml <- time_mlbd <- Kcor_org <- Kcor_ml <- Kcor_mlbd <- AE <- rep(NA, nrep)
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
      time_org[i] <- median(microbenchmark::microbenchmark(Kcor_org[i] <-
      latentcor(X = cbind(x1, x2), types = c("bin", "bin"), method = "original")$R[1, 2], times = 5)$time) / 10^6
      time_ml[i] <- median(microbenchmark::microbenchmark(Kcor_ml[i] <-
      latentcor(X = cbind(x1, x2), types = c("bin", "bin"), ratio = 1)$R[1, 2], times = 5)$time) / 10^6
      time_mlbd[i] <- median(microbenchmark::microbenchmark(Kcor_mlbd[i] <-
      latentcor(X = cbind(x1, x2), types = c("bin", "bin"), ratio = 0.9)$R[1, 2], times= 5)$time) /10^6
    }
    AE <- abs(cbind(Kcor_org - latentRseq[trueR], Kcor_ml - latentRseq[trueR], Kcor_mlbd - latentRseq[trueR],
                    Kcor_ml - Kcor_org, Kcor_mlbd - Kcor_org))
    BB_eval <- c(median(time_org), median(time_ml), median(time_mlbd), colMeans(AE), apply(AE, 2, max))
  }

BB_eval_3d <- array(NA, c(length(latentRseq), length(zratioseq), 13))
for (j in 1:length(latentRseq)) {
  BB_eval_3d[j, , ] = BB_eval[[j]]
}
save(BB_eval_3d, file = "BB_eval.rda", compress = "xz")
