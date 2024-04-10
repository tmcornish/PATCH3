#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --account=m1727

cd ..
srun python -u split_metadata.py
cd slurm_scripts