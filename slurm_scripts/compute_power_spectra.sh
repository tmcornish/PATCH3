#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH --time=24:00:00
#SBATCH -J computePowerSpectra
#SBATCH --mail-user=thomas.cornish@physics.ox.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -A m1727

cd ..

export OMP_NUM_THREADS=256
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

srun -n 1 -c 256 --cpu_bind=cores python -u compute_power_spectra.py

cd slurm_scripts