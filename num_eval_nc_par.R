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
zratioseq <- seq(0.1, 0.8, by = 0.1)
##### check NC
type1 <- "ternary"; type2 <- "continuous"
# the computation results will be saved

registerDoFuture()
plan(multicore, workers = 80)

NC_eval <-
foreach (trueR = 1:length(latentRseq)) %:%
  foreach (zrate = 1:length(zratioseq), .combine = rbind) %dopar% {
    zrate2 <- zratioseq[zrate] + (1 - zratioseq[zrate]) / 2
    # initialize for every combination
    time_org <- time_ml <- time_mlbd <- Kcor_org <- Kcor_ml <- Kcor_mlbd <- AE <- rep(NA, nrep)
    for(i in 1:nrep){
      # generate bivariate normal
      z <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, latentRseq[trueR], latentRseq[trueR], 1), nrow=2))
      z1shift <- quantile(z[, 1], c(zratioseq[zrate], zrate2))
      z1 <- z[, 1] - z1shift[1]
      z2 <- z[, 2] - z1shift[1]
      u1 <- ifelse(z1 > z1shift[2], 2, 1)
      u1[z1 <= 0] = 0          
      u2 <- z2
      # didn't apply any transformation.
      x1 <- u1
      x2 <- u2
      time_org[i] <- median(microbenchmark::microbenchmark(Kcor_org[i] <-
      latentcor(X = cbind(x1, x2), types = c("ter", "con"), method = "original")$R[1, 2], times = 5)$time) / 10^6
      time_ml[i] <- median(microbenchmark::microbenchmark(Kcor_ml[i] <-
      latentcor(X = cbind(x1, x2), types = c("ter", "con"), ratio = 1)$R[1, 2], times = 5)$time) / 10^6
      time_mlbd[i] <- median(microbenchmark::microbenchmark(Kcor_mlbd[i] <-
      latentcor(X = cbind(x1, x2), types = c("ter", "con"), ratio = 0.9)$R[1, 2], times= 5)$time) / 10^6
    }
    AE <- abs(cbind(Kcor_org - latentRseq[trueR], Kcor_ml - latentRseq[trueR], Kcor_mlbd - latentRseq[trueR],
                    Kcor_ml - Kcor_org, Kcor_mlbd - Kcor_org))
    NC_eval <- c(median(time_org), median(time_ml), median(time_mlbd), colMeans(AE), apply(AE, 2, max))
  }

NC_eval_3d <- array(NA, c(length(latentRseq), length(zratioseq), 13))
for (j in 1:length(latentRseq)) {
  NC_eval_3d[j, , ] = NC_eval[[j]]
}
save(NC_eval_3d, file = "NC_eval.rda", compress = "xz")
