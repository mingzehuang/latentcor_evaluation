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

# For TT Case
TTvalue = function (tau_grid, d1_grid, d2_grid) {
  l_tau_grid = length(tau_grid); l_d1_grid = length(d1_grid); l_d2_grid = length(d2_grid)
  TTvalue = array(NA, c(l_tau_grid, l_d1_grid, l_d2_grid))
  registerDoFuture()
  plan(multicore, workers = 80)
  value_list <-
    foreach (i = 1:l_tau_grid) %:%
      foreach (j = 1:l_d1_grid, .combine = rbind) %dopar% {
        value = rep(NA, l_d2_grid)
        for (k in 1:l_d2_grid) {
          zratio1 = d1_grid[j]; zratio2 = d2_grid[k]
          tau = tau_grid[i] * bound_switch(comb = "22", zratio1 = zratio1, zratio2 = zratio2)
          value[k] = r_sol(K = tau, zratio1 = zratio1, zratio2 = zratio2, comb = "22", tol = 1e-8, ratio = 0)
        }
        value_list <- value
      }
  for (i in 1:l_tau_grid) {
    TTvalue[i, , ] = matrix(as.integer(10^7 * value_list[[i]]), l_d1_grid, l_d2_grid)
  }
  return (TTvalue)
}

# grid values that used to create precomputed values.
d1 <- d2 <- log10(seq(1, 10^0.99, length = 50))
tau <- seq(-0.99, 0.99, by = 0.01) # "by" increased from 0.005 to 0.01.

# create grid input for ipol
TTipolgrid <- list(tau, d1, d2)
# interpolation.
TTipol <- chebpol::ipol(TTvalue(tau_grid = tau, d1_grid = d1, d2_grid = d2), grid = TTipolgrid, method = "multilin")
save(TTipol, file = "TT_grid_mixedCCA.rda", compress = "xz")

# grid values that used to create precomputed values
tau_grid <- round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6) * 2 - 1
d1_grid <- d2_grid <- round(pnorm(seq(.12, 1.2, by =.03), sd = .5), 6) * 2 - 1

# create grid input for ipol
TTipolgrid <- list(tau_grid, d1_grid, d2_grid)
# interpolation.
TTipol <- chebpol::ipol(TTvalue(tau_grid = tau_grid, d1_grid = d1_grid, d2_grid = d2_grid), grid = TTipolgrid, method = "multilin")
save(TTipol, file = "TT_grid.rda", compress = "xz")
