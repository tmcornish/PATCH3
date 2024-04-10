#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --account=m1727

cd ../data_query
srun python -u get_data.py
cd ../slurm_scripts