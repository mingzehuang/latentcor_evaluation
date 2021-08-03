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

# For NC Case
NCvalue = function (tau_grid, d11_grid, d12_grid) {
  l_tau_grid = length(tau_grid); l_d11_grid = length(d11_grid); l_d12_grid = length(d12_grid)
  NCvalue = array(NA, c(l_tau_grid, l_d11_grid, l_d12_grid))
  registerDoFuture()
  plan(multicore, workers = 80)
  value_list <-
    foreach (i = 1:l_tau_grid) %:%
      foreach (j = 1:l_d11_grid, .combine = rbind) %dopar% {
        value <- rep(NA, l_d12_grid)
        for (k in 1:l_d12_grid) {
          zratio1 = c(d11_grid[j] * d12_grid[k], d12_grid[k]); zratio2 = NA
          tau = tau_grid[i] * bound_switch(comb = "30", zratio1 = zratio1, zratio2 = zratio2)
          value[k] = r_sol(K = tau, zratio1 = zratio1, zratio2 = zratio2, comb = "30", tol = 1e-8, ratio = 0)
        }
        value_list <- value
      }
  for (i in 1:l_tau_grid) {
    NCvalue[i, , ] = matrix(as.integer(round(10^7 * value_list[[i]], 0)), l_d11_grid, l_d12_grid)
  }
  return (NCvalue)
}

# grid values that used to create precomputed values
tau_grid <- round(pnorm(seq(-1.2, 1.2, by =.06), sd = .5), 6) * 2 - 1
d11_grid <- d12_grid <- round(pnorm(seq(-2.1, 2.1, by =.3)), 6)

# create grid input for ipol
NCipolgrid <- list(tau_grid, d11_grid, d12_grid)
# interpolation.
NCipol <- chebpol::ipol(NCvalue(tau_grid = tau_grid, d11_grid = d11_grid, d12_grid = d12_grid), grid = NCipolgrid, method = "multilin")
save(NCipol, file = "NC_grid.rda", compress = "xz")
