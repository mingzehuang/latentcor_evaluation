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

#for NN Case
# grid values that used to create precomputed values
tauhat_grid <- seq(-1, 1, by = 0.1) # "by" increased from 0.005 to 0.01.
d11_grid <- d12_grid <- d21_grid <- d22_grid <- seq(0.1, 0.9, by = 0.1)
l_tauhat_grid <- length(tauhat_grid); l_d11_grid <- length(d11_grid); l_d12_grid <- length(d12_grid)
l_d21_grid <- length(d21_grid); l_d22_grid <- length(d22_grid)
NNvalue <- array(NA, c(l_tauhat_grid, l_d11_grid, l_d12_grid, l_d21_grid, l_d22_grid))

cl <- makeForkCluster(detectCores(logical=T))
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tauhat_grid) %:%
    foreach (j = 1:l_d11_grid) %dopar% {
      value = array(NA, c(l_d12_grid, l_d21_grid, l_d22_grid))
      for (k in 1:l_d12_grid) {
        for (l in 1:l_d21_grid) {
          for (m in 1:l_d22_grid) {
            zratio1 = c(d11_grid[j] * d12_grid[k], d12_grid[k]);
            zratio2 = c(d21_grid[l] * d22_grid[m], d22_grid[m])
          f1 <- function(r)(bridgeF_nn(r, zratio1 = zratio1, zratio2 = zratio2)
                               - tauhat_grid[i] * bound_nn(zratio1 = zratio1, zratio2 = zratio2))^2
            op <- tryCatch(optimize(f1, lower = -0.999, upper = 0.999, tol = 1e-3)[1], error = function(e) 100)
            if(op == 100) {
              warning("Optimize returned error one of the pairwise correlations, returning NA")
              value[k, l, m] <- NA
            } else {
              value[k, l, m] <- unlist(op)
            }
          }
        }
      }
      value_list <- value
    }
stopCluster(cl)

for (i in 1:l_tauhat_grid) {
  for (j in 1:l_d11_grid) {
    NNvalue[i, j, , , ] = array(as.integer(round(10^7 * value_list[[i]][[j]], 0)), c(l_d12_grid, l_d21_grid, l_d22_grid))
  }
}

# create grid input for ipol
NNipolgrid <- list(tauhat_grid, d11_grid, d12_grid, d21_grid, d22_grid)
# interpolation.
NNipol <- chebpol::ipol(NNvalue, grid = NNipolgrid, method = "multilin")

save(NNipol, file = "NN_grid.rda", compress = "xz")
