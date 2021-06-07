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

# For TB Case
TBvalue = function (tau_grid, d1_grid, d2_grid) {
  l_tau_grid = length(tau_grid); l_d1_grid = length(d1_grid); l_d2_grid = length(d2_grid)
  TBvalue = array(NA, c(l_tau_grid, l_d1_grid, l_d2_grid))
  registerDoFuture()
  plan(multicore, workers = 80)
  value_list <-
    foreach (i = 1:l_tau_grid) %:%
      foreach (j = 1:l_d1_grid, .combine = rbind) %dopar% {
        value = rep(NA, l_d2_grid)
        for (k in 1:l_d2_grid) {
          zratio1 = matrix(d1_grid[j], nrow = 1); zratio2 = matrix(d2_grid[k], nrow = 1)
          tau = tau_grid[i] * bound_tb(zratio1 = zratio1, zratio2 = zratio2)
          value[k] = r_sol(type1 = "trunc", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 1e-6)
        }
        value_list <- value
      }
  for (i in 1:l_tau_grid) {
    TBvalue[i, , ] = matrix(as.integer(10^7 * value_list[[i]]), l_d1_grid, l_d2_grid)
  }
  return (TBvalue)
}

# grid values that used to create precomputed values
d1 <- log10(seq(1, 10^0.99, length = 50))
d2 <- seq(0.01, 0.99, length.out = 50)
tau1 <- c(seq(-0.5, -0.1, by = 0.007), seq(-0.095, -0.001, by = 0.005))
tau <- c(tau1, 0, rev(-tau1))

# create grid input for ipol
TBipolgrid <- list(tau, d1, d2)
# interpolation.
TBipol <- chebpol::ipol(TBvalue(tau_grid = tau, d1_grid = d1, d2_grid = d2), grid = TBipolgrid, method = "multilin")
save(TBipol, file = "TB_grid_mixedCCA.rda", compress = "xz")

# grid values that used to create precomputed values
tau_grid <- round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6) * 2 - 1
d1_grid <- log(seq(1.1, 1000^0.99, length = 20),1000)
d2_grid <- round(pnorm(seq(-1.2, 1.2, by =.12), sd = .5), 6)

# create grid input for ipol
TBipolgrid <- list(tau_grid, d1_grid, d2_grid)
# interpolation.
TBipol <- chebpol::ipol(TBvalue(tau_grid = tau_grid, d1_grid = d1_grid, d2_grid = d2_grid), grid = TBipolgrid, method = "multilin")
save(TBipol, file = "TB_grid.rda", compress = "xz")

