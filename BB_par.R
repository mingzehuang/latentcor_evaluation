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

# For BB Case
# grid values that used to create precomputed values
tau_grid <- round(pnorm(seq(-2.1, 2.1, by = .15)), 6) * 2 - 1 
d1_grid <- d2_grid <- round(pnorm(seq(-2.1, 2.1, by =.3)), 6)
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid); l_d2_grid <- length(d2_grid)
BBvalue <- array(NA, c(l_tau_grid, l_d1_grid, l_d2_grid))

registerDoFuture()
plan(multicore, workers = 80)
value_list <-
  foreach (i = 1:l_tau_grid) %:%
    foreach (j = 1:l_d1_grid, .combine = rbind) %dopar% {
      value <- rep(NA, l_d2_grid)
      for (k in 1:l_d2_grid) {
        zratio1 = matrix(d1_grid[j], nrow = 1); zratio2 = matrix(d2_grid[k], nrow = 1)
        tau = tau_grid[i] * bound_bb(zratio1, zratio2)
        value[k] = r_sol(type1 = "binary", type2 = "binary", tau = tau, zratio1 = zratio1, zratio2 = zratio2, tol = 1e-6)
      }
      value_list <- value
    }

for (i in 1:l_tau_grid) {
      BBvalue[i, , ] = matrix(as.integer(10^7 * value_list[[i]]), l_d1_grid, l_d2_grid)
}

# create grid input for ipol
BBipolgrid <- list(tau_grid, d1_grid, d2_grid)

# interpolation.
BBipol <- chebpol::ipol(BBvalue, grid = BBipolgrid, method = "multilin")
save(BBipol, file = "BB_grid.rda", compress = "xz")


