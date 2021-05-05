load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/BC_grid.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/BB_grid.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/TC_grid.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/TB_grid.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/TT_grid.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/NC_grid.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/NB_grid.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/NN_grid.rda")

save(BCipol, BBipol, TCipol, TBipol, TTipol, NCipol, NBipol, NNipol, file = "all_grid.rda",
     compress = "xz")
