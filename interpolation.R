#Library packages you need for internal functions
library(foreach)
library(doFuture)
library(stats)
library(pcaPP)
library(fMultivar)
library(mnormt)
library(chebpol)
library(Matrix)
library(utils)
library(MASS)
source("/scratch/user/sharkmanhmz/latentcor_git/latentcor/R/internal.R")

value = function(evalfun, grid_list, cores = 80) {
  grid_all = expand.grid(grid_list)
  registerDoFuture()
  plan(multicore, workers = cores)
  value_vector =
    foreach (j = 1:nrow(grid_all), .combine = c) %dopar% {
      grid_input = grid_all[j, ]
      value_list = evalfun(grid_input = grid_input)
    }
  d_grid = length(grid_list)
  dim_value = NULL
  for (i in 1:d) {
    dim_value = c(dim_value, length(grid_list[[i]]))
  }
  return (array(as.integer(10^7 * value_list), dim = dim_value))
}

evalfun_BC = function(grid_input){
  comb = "10"; zratio1 = grid_input[2]; zratio2 = NULL
  K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
  r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
}

# evalfun_BB = function(grid_input){
#   comb = "11"; zratio1 = grid_input[2]; zratio2 = grid_input[3]
#   K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
#   r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
# }
# 
# evalfun_TC = function(grid_input){
#   comb = "20"; zratio1 = grid_input[2]; zratio2 = NULL
#   K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
#   r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
# }
# 
# evalfun_TB = function(grid_input){
#   comb = "21"; zratio1 = grid_input[2]; zratio2 = grid_input[3]
#   K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
#   r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
# }
# 
# evalfun_TT = function(grid_input){
#   comb = "22"; zratio1 = grid_input[2]; zratio2 = grid_input[3]
#   K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
#   r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
# }
# 
# evalfun_NC = function(grid_input){
#   comb = "30"; zratio1 = c(grid_input[2] * grid_input[3], grid_input[3]); zratio2 = NULL
#   K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
#   r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
# }
# 
# evalfun_NB = function(grid_input){
#   comb = "31"; zratio1 = c(grid_input[2] * grid_input[3], grid_input[3]); zratio2 = grid_input[4]
#   K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
#   r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
# }
# 
# evalfun_NT = function(grid_input){
#   comb = "32"; zratio1 = c(grid_input[2] * grid_input[3], grid_input[3]); zratio2 = grid_input[4]
#   K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
#   r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
# }
# 
# evalfun_NN = function(grid_input){
#   comb = "33"; zratio1 = c(grid_input[2] * grid_input[3], grid_input[3]); zratio2 = grid_input[4:5]
#   K = grid_input[1] * bound_switch(comb = comb, zratio1 = zratio1, zratio2 = zratio2)
#   r_switch(method = "original", K = K, zratio1 = zratio1, zratio2 = zratio2, comb = comb, tol = 1e-8, ratio = 0.9)
# }

BCgrid_list = list(seq(-0.99, 0.99, by = 0.02), seq(0.02, 0.98, by = 0.02))
BCipol = chebpol::ipol(value(evalfun = evalfun_BC, grid_list = BCgrid_list), grid = BCgrid_list, method = "multilin")
save(BCipol, file = "interpolation_BC.rda", compress = "xz")
# BBgrid_list = list(seq(-0.99, 0.99, by = 0.03), seq(0.02, 0.98, by = 0.03), seq(0.02, 0.98, by = 0.03))
# BBipol = chebpol::ipol(value(evalfun = evalfun_BB, grid_list = BBgrid_list), grid = BBgrid_list, method = "multilin")
# 
# TCgrid_list = list(seq(-0.99, 0.99, by = 0.02), seq(0.02, 0.98, by = 0.02))
# TCipol = chebpol::ipol(value(evalfun = evalfun_TC, grid_list = TCgrid_list), grid = TCgrid_list, method = "multilin")
# 
# TBgrid_list = list(seq(-0.99, 0.99, by = 0.03), seq(0.02, 0.98, by = 0.03), seq(0.02, 0.98, by = 0.03))
# TBipol = chebpol::ipol(value(evalfun = evalfun_TB, grid_list = TBgrid_list), grid = TBgrid_list, method = "multilin")
# 
# TTgrid_list = list(seq(-0.99, 0.99, by = 0.03), seq(0.02, 0.98, by = 0.03), seq(0.02, 0.98, by = 0.03))
# TTipol = chebpol::ipol(value(evalfun = evalfun_TT, grid_list = TTgrid_list), grid = TTgrid_list, method = "multilin")
# 
# NCgrid_list = list(seq(-0.99, 0.99, by = 0.03), seq(0.02, 0.98, by = 0.03), seq(0.02, 0.98, by = 0.03))
# NCipol = chebpol::ipol(value(evalfun = evalfun_NC, grid_list = NCgrid_list), grid = NCgrid_list, method = "multilin")
# 
# NBgrid_list = list(round(seq(-0.95, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)))
# NBipol = chebpol::ipol(value(evalfun = evalfun_NB, grid_list = NBgrid_list), grid = NBgrid_list, method = "multilin")
# 
# NTgrid_list = list(round(seq(-0.95, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)), round(seq(0.05, 0.95, by = 0.05)))
# NTipol = chebpol::ipol(value(evalfun = evalfun_NT, grid_list = NTgrid_list), grid = NTgrid_list, method = "multilin")
# 
# NNgrid_list = list(seq(-0.9, 0.9, by = 0.1), seq(0.1, 0.9, by = 0.1), seq(0.1, 0.9, by = 0.1), seq(0.1, 0.9, by = 0.1), seq(0.1, 0.9, by = 0.1))
# NNipol = chebpol::ipol(value(evalfun = evalfun_NN, grid_list = NNgrid_list), grid = NNgrid_list, method = "multilin")
# 
# save(BCipol, BBipol, TCipol, TBipol, TTipol, NCipol, NBipol, NTipol, NNipol, file = "interpolation.rda", compress = "xz")


