#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=256
#SBATCH --time=96:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=normal
#SBATCH --account=m1727

cd ..

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=256

srun python -u compute_power_spectra.py

cd slurm_scripts