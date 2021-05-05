#!/bin/bash
##ENVIRONMENT SETTINGS; CHANGE WITH CAUTION
#SBATCH --export=NONE                #Do not propagate environment
#SBATCH --get-user-env=L             #Replicate login environment
#SBATCH --partition=bigmem

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=NN_par            #Set the job name to "JobExample3"
#SBATCH --time=1:00:00              #Set the wall clock limit to 1 Day and 12hr
#SBATCH --nodes=1                #Request 8 tasks
#SBATCH --ntasks-per-node=80         #Request 2 tasks/cores per node
#SBATCH --mem=800GB                    #Request 4096MB (4GB) per node
#SBATCH --output=NNOut.%j      #Send stdout/err to "Example3Out.[jobID]"

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=sharkmanhmz@tamu.edu    #Send all emails to email_address

#First Executable Line
module load GCC/9.3.0  OpenMPI/4.0.3 R/4.0.0
R CMD BATCH --no-save --no-restore --slave NN_par.R

