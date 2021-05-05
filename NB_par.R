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

# For NB Case
# grid values that used to create precomputed values
tauhat_grid <- seq(-1, 1, by = 0.1) # "by" increased from 0.005 to 0.01.
d11_grid <- d12_grid <- d2_grid <- seq(0.1, 0.9, by = 0.1)
l_tauhat_grid <- length(tauhat_grid); l_d11_grid <- length(d11_grid)
l_d12_grid <- length(d12_grid); l_d2_grid <- length(d2_grid)
NBvalue <- array(NA, c(l_tauhat_grid, l_d11_grid, l_d12_grid, l_d2_grid))

cl <- makeForkCluster(detectCores(logical=T))
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tauhat_grid) %:%
  foreach (j = 1:l_d11_grid) %dopar% {
    value = matrix(NA, nrow = l_d12_grid, ncol = l_d2_grid)
    for (k in 1:l_d12_grid) {
      for (l in 1:l_d2_grid) {
        zratio1 = c(d11_grid[j] * d12_grid[k], d12_grid[k]); zratio2 = d2_grid[l]
        f1 <- function(r)(bridgeF_nb(r, zratio1 = zratio1, zratio2 = zratio2)
                         - tauhat_grid[i] * bound_nb(zratio1 = zratio1, zratio2 = zratio2))^2
        op <- tryCatch(optimize(f1, lower = -0.999, upper = 0.999, tol = 1e-3)[1], error = function(e) 100)
        if(op == 100) {
          warning("Optimize returned error one of the pairwise correlations, returning NA")
          value[k, l] <- NA
        } else {
          value[k, l] <- unlist(op)
        }
      }
    }
    value_list <- value
  }
stopCluster(cl)

for (i in 1:l_tauhat_grid) {
  for (j in 1:l_d11_grid) {
        NBvalue[i, j, , ] = matrix(as.integer(round(10^7 * value_list[[i]][[j]], 0)), l_d12_grid, l_d2_grid)
  }
}

# create grid input for ipol
NBipolgrid <- list(tauhat_grid, d11_grid, d12_grid, d2_grid)
# interpolation.
NBipol <- chebpol::ipol(NBvalue, grid = NBipolgrid, method = "multilin")
save(NBipol, file = "NB_grid.rda", compress = "xz")
