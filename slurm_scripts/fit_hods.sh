#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH --time=48:00:00
#SBATCH -J fitHods
#SBATCH --mail-user=thomas.cornish@physics.ox.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -A m1727

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

cd ..

srun -n 1 -c 256 --cpu_bind=none python -u fit_hods.py

cd slurm_scripts