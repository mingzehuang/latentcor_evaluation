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

# For BC Case
BCvalue = function (tau_grid, d1_grid) {
  l_tau_grid = length(tau_grid); l_d1_grid = length(d1_grid)
  registerDoFuture()
  plan(multicore, workers = 80)
  value_list <-
    foreach (i = 1:l_tau_grid, .combine = rbind) %:%
      foreach (j = 1:l_d1_grid, .combine = c) %dopar% {
        zratio1 = d1_grid[j]; zratio2 = NA
        tau = tau_grid[i] * bound_switch(comb = "10", zratio1 = zratio1, zratio2 = zratio2)
        value_list = r_sol(K = tau, zratio1 = zratio1, zratio2 = zratio2, comb = "10", tol = 1e-8, ratio = 0)
      }
  return (matrix(as.integer(10^7 * value_list), l_tau_grid, l_d1_grid))
}

d1 <- seq(0.01, 0.99, length.out = 50)
tau1 <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
tau <- c(tau1, 0, rev(-tau1))

# create grid input for ipol
BCipolgrid <- list(tau, d1)
# interpolation.
BCipol <- chebpol::ipol(BCvalue(tau_grid = tau, d1_grid = d1), grid = BCipolgrid, method = "multilin")
save(BCipol, file = "BC_grid_mixedCCA.rda", compress = "xz")

# grid values that used to create precomputed values
tau_grid <- round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6) * 2 - 1
d1_grid <- round(pnorm(seq(-1.2, 1.2, by =.12), sd = .5), 6)

# create grid input for ipol
BCipolgrid <- list(tau_grid, d1_grid)
# interpolation.
BCipol <- chebpol::ipol(BCvalue(tau_grid = tau_grid, d1_grid = d1_grid), grid = BCipolgrid, method = "multilin")
save(BCipol, file = "BC_grid.rda", compress = "xz")
