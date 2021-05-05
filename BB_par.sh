#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment
#SBATCH --partition=bigmem

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=BB_par            #Set the job name to
#SBATCH --time=0:10:00               #Set the wall clock limit to 15min
#SBATCH --nodes=1                    #Request 1 node
#SBATCH --ntasks-per-node=80         #Request 8 tasks/cores per node
#SBATCH --mem=40GB                     #Request 8GB per node
#SBATCH --output=BBOut.%j            #Send stdout/err to "Example2Out.[jobID]"

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL                     #Send email on all job events
#SBATCH --mail-user=sharkmanhmz@tamu.edu    #Send all emails

#First Executable Line
module load GCC/9.3.0  OpenMPI/4.0.3 R/4.0.0
R CMD BATCH --no-save --no-restore --slave BB_par.R

