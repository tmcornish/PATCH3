#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=256
#SBATCH --time=24:00:00
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --account=m1727

cd ..

#see if pipeline configured to split metadata
if $(python -c "import config as cf; print(cf.makeMapsFromMetadata.split_by_band)")="True"
then
    #get the band and run the script for each one
    for b in $(python -c "import config as cf; print(' '.join(cf.cf_global.bands))")
    do
        srun --cpu-bind=none python -u make_maps_from_metadata.py $b
    done
else
    #get the list of all bands and run them simultaneously
    b=$(python -c "import config as cf; print(','.join(cf.cf_global.bands))")
    srun --cpu-bind=none python -u make_maps_from_metadata.py $b
fi

cd slurm_scripts