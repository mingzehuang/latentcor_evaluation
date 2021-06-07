#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment
#SBATCH --partition=bigmem
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=all_num_eval            #Set the job name to
#SBATCH --time=01:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --nodes=2                    #Request 1 node
#SBATCH --ntasks-per-node=80         #Request 8 tasks/cores per node
#SBATCH --mem=320GB                     #Request 8GB per node
#SBATCH --output=all_num_eval_Out.%j            #Send stdout/err to "Example2Out.[jobID]"

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL                     #Send email on all job events
#SBATCH --mail-user=sharkmanhmz@tamu.edu    #Send all emails

#First Executable Line
module load GCC/10.2.0 OpenMPI/4.0.5 R/4.0.3
R CMD BATCH --no-save --no-restore --slave num_eval_bc_par.R &
R CMD BATCH --no-save --no-restore --slave num_eval_bb_par.R &
R CMD BATCH --no-save --no-restore --slave num_eval_tc_par.R &
R CMD BATCH --no-save --no-restore --slave num_eval_tb_par.R &
R CMD BATCH --no-save --no-restore --slave num_eval_tt_par.R &
R CMD BATCH --no-save --no-restore --slave num_eval_nc_par.R &
R CMD BATCH --no-save --no-restore --slave num_eval_nb_par.R &
R CMD BATCH --no-save --no-restore --slave num_eval_nt_par.R &
R CMD BATCH --no-save --no-restore --slave num_eval_nn_par.R &
wait
R CMD BATCH --no-save --no-restore --slave all_num_eval.R
wait
rm BC_eval.rda BB_eval.rda TC_eval.rda TB_eval.rda TT_eval.rda NC_eval.rda NB_eval.rda NT_eval.rda NN_eval.rda

