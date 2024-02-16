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
# Expects two arguments plus an optional third:
#   $1: Arguments and their values passed to addqueue.
#   $2: Name of the script being run.
#   $3: (Optional) Arguments to pass to the python script itself.
function submit_job () {
    if [ -f $jobfile ]
    then
        jobID=$(getID)
        addqueue $1 --runafter $jobID $PYEX -u $2 $3 > $jobfile
    else
        addqueue $1 $PYEX -u $2 $3 > $jobfile
    fi    
}


# Function for specifying the conditions of running make_maps_from_metadata.py.
function metamaps_job () {
    #see if pipeline configured to split metadata
    if $($PYEX -c "import config as cf; print(cf.makeMapsFromMetadata.split_by_band)")="True"
    then
        #get the band and run the script for each one
        for b in $($PYEX -c "import config as cf; print(' '.join(cf.cf_global.bands))")
        do
            submit_job $1 $2 $b
        done
    else
        #get the list of all bands and run them simultaneously
        b=$($PYEX -c "import config as cf; print(','.join(cf.cf_global.bands))")
        submit_job $1 $2 $b
    fi
}

##### Uncomment all steps below that you wish to run. #####

##### If jobfile exists from previous run, delete it #####
if [ -f $jobfile ]
then
    rm -f $jobfile
fi

### downloading data
#cd data_query/ && submit_job "-q cmb -m 10" get_data.py; cd ..

### splitting the metadata according to field (and possibly filter)
#submit_job "-q cmb -m 20" split_metadata.py

### applying various cuts to clean the catalogues
#submit_job "-q cmb -m 40" clean_catalogues.py

### making maps from the catalogue data
submit_job "-q cmb -m 40" make_maps_from_catalogue.py

### making maps from the frame metadata
# metamaps_job "-q cmb -n 1x20 -m 5 -s" make_maps_from_metadata.py

### making galaxy count and overdensity maps in tomographic bins
#submit_job "-q cmb -m 40" make_galaxy_maps.py

### computing power spectra
#submit_job "-q cmb -n 1x20 -m 2 compute_power_spectra.py" 

