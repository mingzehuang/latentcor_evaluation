library(stats)
library(MASS)
library(pcaPP)
library(Matrix)
library(fMultivar)
library(mnormt)
library(irlba)
library(chebpol)
library(doParallel)
library(foreach)

source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/bridge.R")

#For TC Case
# grid values that used to create precomputed values
tau_grid <- seq(-1, 1, by = 0.1) # "by" increased from 0.005 to 0.01.
d1_grid <- seq(0.1, 0.9, by = 0.1)
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid)
TCvalue <- matrix(NA, l_tau_grid, l_d1_grid)

cl <- makeForkCluster(detectCores(logical=T)) # Use makeForkCluster() on Linux-based system.
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid, .combine = rbind) %:%
    foreach (j = 1:l_d1_grid, .combine = c) %dopar% {
      zratio1 = d1_grid[j]
      f1 <- function(r)(bridgeF_tc(r, zratio1 = zratio1) - tau_grid[i] * bound_tc(zratio1 = zratio1))^2
      op <- tryCatch(optimize(f1, lower = -0.999, upper = 0.999, tol = 1e-3)[1], error = function(e) 100)
      if(op == 100) {
        warning("Optimize returned error one of the pairwise correlations, returning NA")
        value_list <- NA
      } else {
        value_list <- unlist(op)
      }
    }
stopCluster(cl)

TCvalue = matrix(as.integer(10^7 * value_list), l_tau_grid, l_d1_grid)

# create grid input for ipol
TCipolgrid <- list(tau_grid, d1_grid)
# interpolation.
TCipol <- chebpol::ipol(TCvalue, grid = TCipolgrid, method = "multilin")
save(TCipol, file = "TC_grid.rda", compress = "xz")

