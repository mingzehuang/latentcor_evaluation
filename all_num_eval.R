load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/BC_eval.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/BB_eval.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/TC_eval.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/TB_eval.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/TT_eval.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/NC_eval.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/NB_eval.rda")
load("/scratch/user/sharkmanhmz/latentcor_evaluation_git/latentcor_evaluation/NN_eval.rda")

save(BC_eval_3d, BB_eval_3d, TC_eval_3d, TB_eval_3d, TT_eval_3d, NC_eval_3d, NB_eval_3d, NN_eval_3d,
     file = "all_eval.rda", compress = "xz")
