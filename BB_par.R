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

# For BB Case
# grid values that used to create precomputed values
tau_grid <- seq(-1, 1, by = 0.1)
d1_grid <- d2_grid <- seq(0.1, 0.9, by = 0.1)
l_tau_grid <- length(tau_grid); l_d1_grid <- length(d1_grid); l_d2_grid <- length(d2_grid)
BBvalue <- array(NA, c(l_tau_grid, l_d1_grid, l_d2_grid))

cl <- makeForkCluster(detectCores(logical=T))
registerDoParallel(cl)
value_list <-
  foreach (i = 1:l_tau_grid) %:%
    foreach (j = 1:l_d1_grid, .combine = rbind) %dopar% {
      value <- rep(NA, l_d2_grid)
      for (k in 1:l_d2_grid) {
        zratio1 = d1_grid[j]; zratio2 = d2_grid[k]
        f1 <- function(r)(bridgeF_bb(r, zratio1 = zratio1, zratio2 = zratio2) - tau_grid[i] * bound_bb(zratio1, zratio2))^2
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

for (i in 1:l_tau_grid) {
      BBvalue[i, , ] = matrix(as.integer(10^7 * value_list[[i]]), l_d1_grid, l_d2_grid)
}

# create grid input for ipol
BBipolgrid <- list(tau_grid, d1_grid, d2_grid)

# interpolation.
BBipol <- chebpol::ipol(BBvalue, grid = BBipolgrid, method = "multilin")
save(BBipol, file = "BB_grid.rda", compress = "xz")


