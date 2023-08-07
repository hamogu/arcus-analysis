#!/bin/bash
#SBATCH --job-name=arf_rmf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-58
python ../RAeff_single_simulation.py /work/submit/hgunther/ARCUS/arfrmf/ ${SLURM_ARRAY_TASK_ID} \
  --wave_lo=1.5 --wave_hi=60. --wave_step=1.0 \
  -n 500000 --channels 1