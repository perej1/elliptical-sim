#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --mem=2G
#SBATCH --output=output.log
#SBATCH --array=1-2

export OMP_NUM_THREADS=1
module load r/4.1.1-python3
i=$SLURM_ARRAY_TASK_ID
arg=$(sed -n "${i} p" sim-args.txt)
srun Rscript ${arg}
