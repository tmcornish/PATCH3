#!/bin/bash

#######################################################################
# Configure and run this file when running the pipeline on glamdring. #
#######################################################################

PIPEDIR=$(pwd)
PYEX=$(which python)    #get the Python executable path
jobfile="$PIPEDIR/prevjob.txt"   #file to which job ID will be output

### If jobfile exists from previous run, delete it ###
if [ -f $jobfile ]
then
    rm -f $jobfile
fi

# Function for retrieving previous job ID from file
function getID () {
    tail -3 $jobfile | head -1 | awk -F'python-' '{ print $NF }' | awk -F'.out' '{ print $1 }'
}


# Function to to submit a bash-scripted job to the queue.
# Expects two arguments plus an optional third:
#   $1: Arguments and their values passed to addqueue.
#   $2: Name of the script being run.
#   $3: (Optional) Arguments to pass to the bash script itself.
function submit_job () {
    if [ -f $jobfile ]
    then
        jobID=$(getID)
        addqueue $1 --runafter $jobID $2 $3 > $jobfile
    else
        addqueue $1 $2 $3 > $jobfile
    fi    
}


# Function to to submit a Python job to the queue.
# Expects two arguments plus an optional third:
#   $1: Arguments and their values passed to addqueue.
#   $2: Name of the script being run.
#   $3: (Optional) Arguments to pass to the python script itself.
function submit_pyjob () {
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
    #see if previous job is running
    if [ -f $jobfile ]
    then
        jobID=$(getID)
        runafter="--runafter $jobID"
    else
        runafter=""
    fi

    #see if pipeline configured to split metadata
    if [ $($PYEX -c "import config as cf; print(cf.makeMapsFromMetadata.split_by_band)")="True" ]
    then
        #get the band and run the script for each one
        for b in $($PYEX -c "import config as cf; print(' '.join(cf.cf_global.bands))")
        do
            addqueue $1 "$runafter" $PYEX -u $2 $b > $jobfile
        done
    else
        #get the list of all bands and run them simultaneously
        b=$($PYEX -c "import config as cf; print(','.join(cf.cf_global.bands))")
        addqueue $1 "$runafter" $PYEX -u $2 $b > $jobfile
    fi
}


# Function for submitting compute_power_spectra.py to the queue.
function power_spectra_job () {
    #see if previous job is running
    if [ -f $jobfile ]
    then
        jobID=$(getID)
        runafter="--runafter $jobID"
    else
        runafter=""
    fi

    for p in $($PYEX -c "import config as cf; print(' '.join(cf.computePowerSpectra.get_bin_pairings()[1]))")
    do
        #submit the job to the queue
        addqueue -q cmb -n 1x$1 -m $2 "$runafter" -s $PYEX $3 $p > $jobfile
    done
}

##### Uncomment all steps below that you wish to run. #####


### downloading data
#cd data_query/ && submit_pyjob "-q cmb -m 10" get_data.py; cd ..

### splitting the metadata according to field (and possibly filter)
#submit_pyjob "-q cmb -m 20" split_metadata.py

### applying various cuts to clean the catalogues
#submit_pyjob "-q cmb -m 40" clean_catalogues.py

### making maps from the catalogue data
#submit_pyjob "-q cmb -m 40" make_maps_from_catalogue.py

### making maps from the frame metadata
#metamaps_job "-q cmb -n 1x24 -m 7 -s" make_maps_from_metadata.py

### making galaxy count and overdensity maps in tomographic bins
#submit_pyjob "-q cmb -m 40" make_galaxy_maps.py

### computing power spectra; function takes as arguments...
###     1: number of cores to use
###     2: memory per CPU
###     3: name of the python script to run
power_spectra_job 24 7 compute_power_spectra.py

