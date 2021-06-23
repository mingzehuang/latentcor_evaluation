library(MASS)
library(microbenchmark)
library(foreach)
library(Matrix)
library(chebpol)
library(pcaPP)
library(doRNG)
library(doFuture)

load("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/sysdata.rda")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/bridge.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/KendallTau.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/fromKtoR.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/estR.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/fromZtoX.R")
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/GenData.R")
# setup for 100 replication.
nrep = 100
# sample size
n = 100
##### check BC
type1 = "binary"; type2 = "continuous"
# will test 19 latent r and 9 zero proportion values.
latentRseq = seq(-0.9, 0.9, by = 0.1)
zratioseq = seq(0.1, 0.9, by = 0.1)
types = c("continuous", "binary", "trunc", "ternary")
methods = c("original", "ml", "approx")
indseq = expand.grid(nrep, latentRseq, zratioseq, types, types, methods)
registerDoFuture()
plan(multicore, workers = 80)
value =
foreach (ind = 1:nrow(indseq)) {
  simdata = GenData(n = n, type1 = type1, type2 = type2, rho = latentRseq[trueR], q1 = zratioseq[zrate])
  X1 = simdata$X1; X2 = simdata$X2
  time = median(microbenchmark::microbenchmark(estimate = estR(X1 = x1, X2 = x2, type1 = type1, type2 = type2, method = methods[m])$R12, times = 5)$time) / 10^6
  value = c(time, estimate)
}
out = array(NA, c(length(latentRseq), length(zratioseq), length(methods), (length(methods)) ^ 3 * (length(methods) - 1) ^ 2))
for (j in 1:nrow(indseq)) {
  out[indseq[j, 1], indseq[j, 2], indseq[j, 3], ] = value[j]
}
save(out, file = "out.rda", compress = "xz")