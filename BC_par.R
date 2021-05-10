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

# For BC Case
# grid values that used to create precomputed values

tau_grid <- round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6) * 2 - 1
d1_grid <- round(pnorm(seq(-1.2, 1.2, by =.12), sd = .5), 6)
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid)
BCvalue <- matrix(NA, l_tau_grid, l_d1_grid)

registerDoFuture()
plan(multicore, workers = 80)

value_list <-
  foreach (i = 1:l_tau_grid, .combine = rbind) %:%
    foreach (j = 1:l_d1_grid, .combine = c) %dopar% {
      zratio1 = matrix(d1_grid[j], nrow = 1)
      tau = tau_grid[i] * bound_bc(zratio1 = zratio1)
      value_list = r_sol(type1 = "binary", type2 = "continuous", tau = tau, zratio1 = zratio1, zratio2 = NULL, tol = 1e-6)
    }

BCvalue = matrix(as.integer(10^7 * value_list), l_tau_grid, l_d1_grid)

# create grid input for ipol
BCipolgrid <- list(tau_grid, d1_grid)
# interpolation.
BCipol <- chebpol::ipol(BCvalue, grid = BCipolgrid, method = "multilin")
save(BCipol, file = "BC_grid.rda", compress = "xz")


