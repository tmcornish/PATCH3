#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:30:00
#SBATCH --constraint=cpu
#SBATCH --qos=debug
#SBATCH --account=m1727

cd ..
srun python -u plot_power_spectra.py
cd slurm_scripts