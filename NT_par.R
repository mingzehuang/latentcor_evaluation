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
# For NT Case
# grid values that used to create precomputed values
tau_grid <- round(pnorm(seq(-1.8, 1.8, by =.15), sd = .8), 6) * 2 - 1
d11_grid <- d12_grid <- round(pnorm(seq(-1.8, 1.8, by =.3), sd = .8), 6)
d2_grid <- round(pnorm(seq(.12, 1.2, by =.06), sd = .5), 6) * 2 - 1
l_tau_grid <- length(tau_grid); l_d11_grid <- length(d11_grid)
l_d12_grid <- length(d12_grid); l_d2_grid <- length(d2_grid)
NTvalue <- array(NA, c(l_tau_grid, l_d11_grid, l_d12_grid, l_d2_grid))

registerDoFuture()
plan(multicore, workers = 80)

value_list <-
  foreach (i = 1:l_tau_grid) %:%
  foreach (j = 1:l_d11_grid) %dopar% {
    value = matrix(NA, nrow = l_d12_grid, ncol = l_d2_grid)
    for (k in 1:l_d12_grid) {
      for (l in 1:l_d2_grid) {
        zratio1 = matrix(c(d11_grid[j] * d12_grid[k], d12_grid[k]), nrow = 1); zratio2 = matrix(d2_grid[l], nrow = 1)
        tau = tau_grid[i] * bound_nt(zratio1 = zratio1, zratio2 = zratio2)
        value[k, l] = r_sol(type1 = "ternary", type2 = "trunc", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 1e-6)  
      }
    }
    value_list <- value
  }

for (i in 1:l_tau_grid) {
  for (j in 1:l_d11_grid) {
        NTvalue[i, j, , ] = matrix(as.integer(round(10^7 * value_list[[i]][[j]], 0)), l_d12_grid, l_d2_grid)
  }
}

# create grid input for ipol
NTipolgrid <- list(tau_grid, d11_grid, d12_grid, d2_grid)
# interpolation.
NTipol <- chebpol::ipol(NTvalue, grid = NTipolgrid, method = "multilin")
save(NTipol, file = "NT_grid.rda", compress = "xz")
