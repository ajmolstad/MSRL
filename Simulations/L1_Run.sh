#!/bin/bash

#SBATCH -o Results/L1_%a.Rout
#SBATCH --array=1-3000
#SBATCH --mail-user=amolstad@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 12:00:00

export OMP_NUM_THREADS=1

module load R/3.6

R CMD BATCH --vanilla Simulations_Main.R  Results/L1_${SLURM_ARRAY_TASK_ID}.Rout
