#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4G
#SBATCH --output=output/out_%A_%a.out
#SBATCH --array=1-108

export OMP_NUM_THREADS=1
module load r/4.1.1-python3
i=$SLURM_ARRAY_TASK_ID
arg=$(sed -n "${i} p" sim-args.txt)
srun Rscript ${arg}
