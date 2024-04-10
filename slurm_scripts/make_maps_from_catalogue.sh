#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=normal
#SBATCH --account=m1727

cd ..
srun python -u make_maps_from_catalogue.py
cd slurm_scripts