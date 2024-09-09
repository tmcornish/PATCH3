#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q debug
#SBATCH --time=00:30:00
#SBATCH -J plotPowerSpectra
#SBATCH --mail-user=thomas.cornish@physics.ox.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -A m1727

cd ..
srun python -u plot_power_spectra.py
cd slurm_scripts