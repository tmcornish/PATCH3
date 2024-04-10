#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=normal
#SBATCH --account=m1727

cd ..
srun python -u clean_catalogues.py
cd slurm_scripts