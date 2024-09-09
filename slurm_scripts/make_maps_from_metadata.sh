#!/bin/bash
#SBATCH -N 1
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH --time=24:00:00
#SBATCH -J makeMapsFromMetadata
#SBATCH --mail-user=thomas.cornish@physics.ox.ac.uk
#SBATCH --mail-type=ALL
#SBATCH -A m1727

cd ..

#see if pipeline configured to split metadata
if [ $(python -c "import config as cf; print(cf.makeMapsFromMetadata.split_by_band)")="True" ]
then
    #get the band and run the script for each one
    for b in $(python -c "import config as cf; print(' '.join(cf.cf_global.bands))")
    do
        srun -n 1 -c 256 --cpu-bind=none python -u make_maps_from_metadata.py $b
    done
else
    #get the list of all bands and run them simultaneously
    b=$(python -c "import config as cf; print(','.join(cf.cf_global.bands))")
    srun -n 1 -c 256 --cpu-bind=none python -u make_maps_from_metadata.py $b
fi

cd slurm_scripts