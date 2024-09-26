#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH --time=00:30:00
#SBATCH -J makeMapsFromCatalogue
#SBATCH --mail-user=thomas.cornish@physics.ox.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -A m1727

cd ..
srun python -u make_maps_from_catalogue.py
cd slurm_scripts