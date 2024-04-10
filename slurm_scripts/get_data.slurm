#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=normal
#SBATCH --account=m1727

PYEX=$(which python)

cd data_query
srun python -u get_data.py
cd ..