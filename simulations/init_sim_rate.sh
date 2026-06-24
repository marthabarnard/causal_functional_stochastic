#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=25
#SBATCH --array=1-240
#SBATCH --mem-per-cpu=10000M
#SBATCH -t 5:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=barna126@umn.edu
#SBATCH --job-name=function_sim_new_rate
#SBATCH -A julianw_queue
#SBATCH -o /users/9/barna126/scalar_sequence2/logs_rate_Q/%A_%a.out
#SBATCH -e /users/9/barna126/scalar_sequence2/logs_rate_Q/%A_%a.err

date
path=/users/9/barna126
cd $path/scalar_sequence2/
module load R/4.4.0-openblas-rocky8
R CMD BATCH --no-save --no-restore $path/scalar_sequence2/code/run_sim_rate.R $path/scalar_sequence2/logs_rate_Q/job_$SLURM_ARRAY_TASK_ID.txt
