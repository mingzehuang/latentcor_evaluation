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
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/bridge.R")

#For TC Case
# grid values that used to create precomputed values
tau_grid <- round(pnorm(seq(-2.1, 2.1, by = .15)), 6) * 2 - 1
d1_grid <- round(pnorm(seq(.15, 2.1, by = .15)), 6) * 2 - 1
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid)
TCvalue <- matrix(NA, l_tau_grid, l_d1_grid)

registerDoFuture()
plan(multicore, workers = 80)
value_list <-
  foreach (i = 1:l_tau_grid, .combine = rbind) %:%
    foreach (j = 1:l_d1_grid, .combine = c) %dopar% {
      zratio1 = matrix(d1_grid[j], nrow = 1)
      tau = tau_grid[i] * bound_tc(zratio1 = zratio1)
      value_list = r_sol(type1 = "trunc", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, tol = 1e-6)
    }

TCvalue = matrix(as.integer(10^7 * value_list), l_tau_grid, l_d1_grid)

# create grid input for ipol
TCipolgrid <- list(tau_grid, d1_grid)
# interpolation.
TCipol <- chebpol::ipol(TCvalue, grid = TCipolgrid, method = "multilin")
save(TCipol, file = "TC_grid.rda", compress = "xz")

