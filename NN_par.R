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

#for NN Case
NNvalue = function (tau_grid, d11_grid, d12_grid, d21_grid, d22_grid) {
  l_tau_grid = length(tau_grid); l_d11_grid = length(d11_grid); l_d12_grid = length(d12_grid)
  l_d21_grid = length(d21_grid); l_d22_grid = length(d22_grid)
  NNvalue = array(NA, c(l_tau_grid, l_d11_grid, l_d12_grid, l_d21_grid, l_d22_grid))
  registerDoFuture()
  plan(multicore, workers = 80)
  value_list <-
    foreach (i = 1:l_tau_grid) %:%
      foreach (j = 1:l_d11_grid) %dopar% {
        value = array(NA, c(l_d12_grid, l_d21_grid, l_d22_grid))
        for (k in 1:l_d12_grid) {
          for (l in 1:l_d21_grid) {
            for (m in 1:l_d22_grid) {
              zratio1 = matrix(c(d11_grid[j] * d12_grid[k], d12_grid[k]), nrow = 1)
              zratio2 = matrix(c(d21_grid[l] * d22_grid[m], d22_grid[m]), nrow = 1)
              tau = tau_grid[i] * bound_nn(zratio1 = zratio1, zratio2 = zratio2)
              value[k, l, m] = r_sol(type1 = "ternary", type2 = "ternary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 1e-6)
            }
          }
        }
        value_list <- value
      }
  for (i in 1:l_tau_grid) {
    for (j in 1:l_d11_grid) {
      NNvalue[i, j, , , ] = array(as.integer(round(10^7 * value_list[[i]][[j]], 0)), c(l_d12_grid, l_d21_grid, l_d22_grid))
    }
  }
  return (NNvalue)
}

# grid values that used to create precomputed values
tau_grid <- round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6) * 2 - 1
d11_grid <- d12_grid <- d21_grid <- d22_grid <- round(pnorm(seq(-1.8, 1.8, by =.3), sd = .8), 6)

# create grid input for ipol
NNipolgrid <- list(tau_grid, d11_grid, d12_grid, d21_grid, d22_grid)
# interpolation.
NNipol <- chebpol::ipol(NNvalue(tau_grid = tau_grid, d11_grid = d11_grid, d12_grid = d12_grid, d21_grid = d21_grid, d22_grid = d22_grid), grid = NNipolgrid, method = "multilin")
save(NNipol, file = "NN_grid.rda", compress = "xz")
