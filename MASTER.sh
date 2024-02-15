#!/bin/bash

#######################################################################
# Configure and run this file when running the pipeline on glamdring. #
#######################################################################

PIPEDIR=$(pwd)
PYEX=$(which python)    #get the Python executable path
jobfile="$PIPEDIR/prevjob.txt"   #file to which job ID will be output

# Function for retrieving previous job ID from file
function getID () {
    sed -n '1p' $jobfile | awk -F'sh-' '{ print $NF }' | awk -F'.out' '{ print $1 }'
}

# Function to to submit a job to the queue.
# Expects two arguments:
#   $1: Arguments and their values passed to addqueue.
#   $2: Name of the script being run.
function submit_job () {
    if [test -f $jobfile]
    then
        jobID=$(getID)
        addqueue $1 --runafter $jobID $PYEX -u $2 > $jobfile
    else
        addqueue $1 $PYEX -u $2 > $jobfile
    fi    
}




##### Uncomment all steps below that you wish to run. #####


### downloading data
#cd data_query/ && submit_job "-q cmb -m 10" get_data.py; cd ..

### splitting the metadata according to field (and possibly filter)
#submit_job "-q cmb -m 20" split_metadata.py

### applying various cuts to clean the catalogues
#submit_job "-q cmb -m 40" clean_catalogues.py

### making maps from the catalogue data
submit_job "-q cmb -m 40" make_maps_from_catalogue.py

### making maps from the frame metadata
# TODO: figure this out

### making galaxy count and overdensity maps in tomographic bins
#submit_job "-q cmb -m 40" make_galaxy_maps.py

### computing power spectra
#submit_job "-q cmb -n 1x20 -m 2 compute_power_spectra.py" 

