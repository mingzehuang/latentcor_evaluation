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

# For BC Case
# grid values that used to create precomputed values

tau_grid <- seq(-1, 1, by = 0.1)
d1_grid <- seq(0.1, 0.9, by = 0.1)
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid)
BCvalue <- matrix(NA, l_tau_grid, l_d1_grid)

cl <- makeForkCluster(detectCores())
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid, .combine = rbind) %:%
    foreach (j = 1:l_d1_grid, .combine = c) %dopar% {
      zratio1 = d1_grid[j]
      f1 <- function(r)(bridgeF_bc(r, zratio1 = zratio1) - tau_grid[i] * bound_bc(zratio1 = zratio1))^2
      op <- tryCatch(optimize(f1, lower = -0.999, upper = 0.999, tol = 1e-3)[1], error = function(e) 100)
      if(op == 100) {
        warning("Optimize returned error one of the pairwise correlations, returning NA")
        value_list <- NA
      } else {
        value_list <- unlist(op)
      }
    }
stopCluster(cl)

BCvalue = matrix(as.integer(10^7 * value_list), l_tau_grid, l_d1_grid)

# create grid input for ipol
BCipolgrid <- list(tau_grid, d1_grid)
# interpolation.
BCipol <- chebpol::ipol(BCvalue, grid = BCipolgrid, method = "multilin")
save(BCipol, file = "BC_grid.rda", compress = "xz")
print(cl)

