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

# For NC Case
# grid values that used to create precomputed values
tauhat_grid <- seq(-1, 1, by = 0.1) # "by" increased from 0.005 to 0.01.
d11_grid <- d12_grid <- seq(0.1, 0.9, by = 0.1)
l_tauhat_grid <- length(tauhat_grid); l_d11_grid <- length(d11_grid); l_d12_grid <- length(d12_grid)
NCvalue <- array(NA, c(l_tauhat_grid, l_d11_grid, l_d12_grid))

cl <- makeForkCluster(detectCores(logical=T))
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tauhat_grid) %:%
    foreach (j = 1:l_d11_grid, .combine = rbind) %dopar% {
      value <- rep(NA, l_d12_grid)
      for (k in 1:l_d12_grid) {
        zratio1 = c(d11_grid[j] * d12_grid[k], d12_grid[k])
        f1 <- function(r)(bridgeF_nc(r, zratio1 = zratio1) - tauhat_grid[i] * bound_nc(zratio1 = zratio1))^2
        op <- tryCatch(optimize(f1, lower = -0.999, upper = 0.999, tol = 1e-3)[1], error = function(e) 100)
        if(op == 100) {
          warning("Optimize returned error one of the pairwise correlations, returning NA")
          value[k] <- NA
        } else {
          value[k] <- unlist(op)
        }
      }
      value_list <- value
    }
stopCluster(cl)

for (i in 1:l_tauhat_grid) {
      NCvalue[i, , ] = matrix(as.integer(round(10^7 * value_list[[i]], 0)), l_d11_grid, l_d12_grid)
}

# create grid input for ipol
NCipolgrid <- list(tauhat_grid, d11_grid, d12_grid)
# interpolation.
NCipol <- chebpol::ipol(NCvalue, grid = NCipolgrid, method = "multilin")
save(NCipol, file = "NC_grid.rda", compress = "xz")


