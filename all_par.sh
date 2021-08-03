#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment
#SBATCH --partition=bigmem
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=all_par            #Set the job name to
#SBATCH --time=01:00:00               #Set the wall clock limit to 6hr and 30min
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=80         #Request 8 tasks/cores per node
#SBATCH --mem=320GB                     #Request 8GB per node
#SBATCH --output=all_Out.%j            #Send stdout/err to "Example2Out.[jobID]"

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL                     #Send email on all job events
#SBATCH --mail-user=sharkmanhmz@tamu.edu    #Send all emails

#First Executable Line
module load GCC/10.2.0 OpenMPI/4.0.5 R/4.0.3
R CMD BATCH --no-save --no-restore --slave BC_par.R &
R CMD BATCH --no-save --no-restore --slave BB_par.R &
R CMD BATCH --no-save --no-restore --slave TC_par.R &
R CMD BATCH --no-save --no-restore --slave TB_par.R &
R CMD BATCH --no-save --no-restore --slave TT_par.R &
R CMD BATCH --no-save --no-restore --slave NC_par.R &
R CMD BATCH --no-save --no-restore --slave NB_par.R &
R CMD BATCH --no-save --no-restore --slave NT_par.R &
R CMD BATCH --no-save --no-restore --slave NN_par.R &
wait
echo "BC_grid_mixedCCA.rda size is $(find "BC_grid_mixedCCA.rda" -printf "%s") Bytes."
echo "BC_grid.rda size is $(find "BC_grid.rda" -printf "%s") Bytes."
echo "BB_grid_mixedCCA.rda size is $(find "BB_grid_mixedCCA.rda" -printf "%s") Bytes."
echo "BB_grid.rda size is $(find "BB_grid.rda" -printf "%s") Bytes."
echo "TC_grid_mixedCCA.rda size is $(find "TC_grid_mixedCCA.rda" -printf "%s") Bytes."
echo "TC_grid.rda size is $(find "TC_grid.rda" -printf "%s") Bytes."
echo "TB_grid_mixedCCA.rda size is $(find "TB_grid_mixedCCA.rda" -printf "%s") Bytes."
echo "TB_grid.rda size is $(find "TB_grid.rda" -printf "%s") Bytes."
echo "TT_grid_mixedCCA.rda size is $(find "TT_grid_mixedCCA.rda" -printf "%s") Bytes."
echo "TT_grid.rda size is $(find "TT_grid.rda" -printf "%s") Bytes."
echo "NC_grid.rda size is $(find "NC_grid.rda" -printf "%s") Bytes."
echo "NB_grid.rda size is $(find "NB_grid.rda" -printf "%s") Bytes."
echo "NT_grid.rda size is $(find "NT_grid.rda" -printf "%s") Bytes."
echo "NN_grid.rda size is $(find "NN_grid.rda" -printf "%s") Bytes."
R CMD BATCH --no-save --no-restore --slave all_par.R
wait
echo "all_grid.rda size is $(find "all_grid.rda" -printf "%s") Bytes."
rm BC_grid_mixedCCA.rda BC_grid.rda BB_grid_mixedCCA.rda BB_grid.rda TC_grid_mixedCCA.rda TC_grid.rda TB_grid_mixedCCA.rda TB_grid.rda TT_grid_mixedCCA.rda TT_grid.rda NC_grid.rda NB_grid.rda NT_grid.rda NN_grid.rda
