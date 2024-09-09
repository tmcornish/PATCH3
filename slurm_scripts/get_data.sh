#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH --time=48:00:00
#SBATCH -J getData
#SBATCH --mail-user=thomas.cornish@physics.ox.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -A m1727

cd ../data_query
srun python -u get_data.py
cd ../slurm_scripts