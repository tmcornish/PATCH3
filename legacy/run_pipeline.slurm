#!/bin/bash
#SBATCH --job-name=pHSC3      # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # number of processes
#SBATCH --mem-per-cpu=40G         # memory per cpu-core (4G per cpu-core is default)
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --partition=cmb			 # queue on which the process is to be run

srun python -u MASTER.py