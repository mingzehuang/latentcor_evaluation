library(stats)
library(MASS)
library(pcaPP)
library(Matrix)
library(fMultivar)
library(mnormt)
library(irlba)
library(chebpol)
library(foreach)
library(doFuture)
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/internal.R")

#For TC Case
TCvalue = function (tau_grid, d1_grid) {
  l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid)
  registerDoFuture()
  plan(multicore, workers = 80)
  value_list <-
    foreach (i = 1:l_tau_grid, .combine = rbind) %:%
      foreach (j = 1:l_d1_grid, .combine = c) %dopar% {
        zratio1 = d1_grid[j]; zratio2 = NA
        tau = tau_grid[i] * bound_switch(comb = "20", zratio1 = zratio1, zratio2 = zratio2)
        value_list = r_sol(K = tau, zratio1 = zratio1, zratio2 = zratio2, comb = "20", tol = 1e-8, ratio = 0)
      }
  return (matrix(as.integer(10^7 * value_list), l_tau_grid, l_d1_grid))
}

# grid values that used to create precomputed values.
d1 <- log10(seq(1, 10^0.99, length = 50))
tau <- seq(-0.99, 0.99, by = 0.01) # "by" increased from 0.005 to 0.01.

# create grid input for ipol
TCipolgrid <- list(tau, d1)
# interpolation.
TCipol <- chebpol::ipol(TCvalue(tau_grid = tau, d1_grid = d1), grid = TCipolgrid, method = "multilin")
save(TCipol, file = "TC_grid_mixedCCA.rda", compress = "xz")

# grid values that used to create precomputed values
tau_grid <- round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6) * 2 - 1
d1_grid <- round(pnorm(seq(.1, 2.5, by = .1)), 6) * 2 - 1

# create grid input for ipol
TCipolgrid <- list(tau_grid, d1_grid)
# interpolation.
TCipol <- chebpol::ipol(TCvalue(tau_grid = tau_grid, d1_grid = d1_grid), grid = TCipolgrid, method = "multilin")
save(TCipol, file = "TC_grid.rda", compress = "xz")

