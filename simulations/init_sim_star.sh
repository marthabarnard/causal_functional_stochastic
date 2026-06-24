#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=25
#SBATCH --array=1-72
#SBATCH --mem-per-cpu=8000M
#SBATCH -t 4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=barna126@umn.edu
#SBATCH --job-name=run_Qstar_sim
#SBATCH -A julianw_queue
#SBATCH -o /users/9/barna126/scalar_sequence2/logs_Qstar/%A_%a.out
#SBATCH -e /users/9/barna126/scalar_sequence2/logs_Qstar/%A_%a.err

date
path=/users/9/barna126
cd $path/scalar_sequence2/
module load R/4.4.0-openblas-rocky8
R CMD BATCH --no-save --no-restore $path/scalar_sequence2/code/run_sim_star.R $path/scalar_sequence2/logs_Qstar/job_$SLURM_ARRAY_TASK_ID.txt
