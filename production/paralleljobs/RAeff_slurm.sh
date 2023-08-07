#!/bin/bash
#SBATCH --job-name=RAeff
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-389
python RAeff_single_simulation.py /work/submit/hgunther/ARCUS/raysRAeff/ ${SLURM_ARRAY_TASK_ID}