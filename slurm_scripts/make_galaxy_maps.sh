#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --constraint=cpu
#SBATCH --qos=normal
#SBATCH --account=m1727

cd ..
srun python -u make_galaxy_maps.py
cd slurm_scripts