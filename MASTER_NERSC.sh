#!/bin/bash

###################################################################
# Configure and run this file when running the pipeline on NERSC. #
###################################################################


#uncomment any stages you want to run
jobs=(
    #'get_data'
    'split_metadata'
    #'clean_catalogues'
    #'make_maps_from_catalogue'
    #'make_maps_from_metadata'
    #'pca_systematics'
    #'make_galaxy_maps'
    #'compute_power_spectra'
    #'plot_power_spectra'
    #'make_txpipe_inputs'
)



module load conda
conda activate phsc3

cd slurm_scripts

#run first job and get job ID
jid=$(sbatch ${jobs[0]}.sh | cut -d ' ' -f4)

for j in ${jobs[@]:1}
    do
        jid=$(sbatch --dependency=afterok:${jid} ${j}.sh | cut -d ' ' -f4)
    done

cd ..