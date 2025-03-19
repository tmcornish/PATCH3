#!/bin/bash

#######################################################################
# Configure and run this file when running the pipeline on glamdring. #
#######################################################################

#first activate the conda environment
source /mnt/users/tcornish/miniconda/etc/profile.d/conda.sh
conda activate patch3

config_file=$1          #config file passed as argument
PIPEDIR=$(pwd)          #current directory
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
        addqueue $1 --runafter $jobID $PYEX -u $2 $config_file $3 > $jobfile
    else
        addqueue $1 $PYEX -u $2 $config_file $3 > $jobfile
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
    pya="from configuration import PipelineConfig as PC; "
    pyb="cf = PC('$config_file', 'makeMapsFromMetadata'); "
    pyc="print(cf.split_by_band)"
    if [[ $($PYEX -c "$pya$pyb$pyc") == "True" ]]
    then
        #get the band and run the script for each one
        pyc="print(' '.join(cf.bands.all))"
        for b in $($PYEX -c "$pya$pyb$pyc")
        do
            addqueue $1 "$runafter" $PYEX -u $2 $config_file $b > $jobfile
        done
    else
        #get the list of all bands and run them simultaneously
        pyc="print(','.join(cf.bands.all))"
        b=$($PYEX -c "$pya$pyb$pyc")
        addqueue $1 "$runafter" $PYEX -u $2 $config_file $b > $jobfile
    fi
}


# Function for adding MPI python jobs to the queue.
# Expects four arguments plus an optional third:
#   $1: Number of nodes required.
#   $2: Number of CPUs required per node.
#   $3: Memory required per CPU (in GB).
#   $4: Name of the script being run.
#   $5: (Optional) Arguments to pass to the python script itself.
function submit_mpi_pyjob () {
    #see if previous job is running
    if [ -f $jobfile ]
    then
        jobID=$(getID)
        runafter="--runafter $jobID"
    else
        runafter=""
    fi
    #lines comprising the command to run the job
    runfile="./run_mpi_job.sh"
    c0=$'#!/bin/bash\n'
    c1="/usr/local/shared/slurm/bin/srun -N $1 --ntasks-per-node=1 --cpus-per-task=$2 "
    c2="--cpu-bind=cores -m cyclic --mpi=pmix $PYEX -u $4 $config_file $5"
    echo "$c0$c1$c2" > $runfile
    chmod u+x $runfile
    addqueue -s -q redwood -n "$1"x"$2" -m $3 "$runafter" $runfile > $jobfile
}

##### Uncomment all steps below that you wish to run. #####


### downloading data
#cd data_query/ && submit_pyjob "-q cmb -m 10" get_data.py; cd ..

### splitting the metadata according to field (and possibly filter)
#submit_pyjob "-q cmb -m 20" split_metadata.py

### applying various cuts to clean the catalogues
#submit_pyjob "-q cmb -m 40" clean_catalogues.py

### selecting galaxy samples for analysis
#submit_pyjob "-q cmb -m 40" sample_selection.py

### making maps from the frame metadata
#metamaps_job "-q cmb -n 1x24 -m 7 -s" make_maps_from_metadata.py

### making maps from the catalogue data
#submit_pyjob "-q cmb -n 1x10 -m 4 -s" make_maps_from_catalogue.py

### making galaxy count and overdensity maps in tomographic bins
#submit_pyjob "-q cmb -m 40" make_galaxy_maps.py

### combining maps from all fields
#submit_pyjob "-q cmb -m 40" combine_fields.py

### compute n(z) distributions using DIR
#submit_pyjob "-q cmb -n 1x24 -m 5 -s" dir_photozs.py

### compute theory predictions of the power spectra
#submit_pyjob "-q cmb -m 40" theory_predictions.py

### computing power spectra; function takes as arguments...
###     1: number of nodes to use
###     2: number of CPUs to use per node
###     3: memory per CPU (in GB)
###     4: name of the python script to run
submit_mpi_pyjob 5 24 4 compute_power_spectra.py

### calculating covariances
#submit_pyjob "-q cmb -n 1x24 -m 7 -s" covariances.py

### calculating covariances
#submit_pyjob "-q cmb -n 1x24 -m 2 -s" make_sacc_files.py

### fitting HOD models to data
#submit_pyjob "-q cmb -n 1x24 -m 7 -s" fit_hods.py