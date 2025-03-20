# Pipeline for Analysis of Tomographic Clustering with HSC DR3
\*\*\* **WORK IN PROGRESS** \*\*\*
PATCH3 is a pipeline for performing Fourier-space galaxy clustering analyses using Data Release 3 (DR3) from the Hyper Suprime-Cam Subaru Strategic Program (HSC-SSP). NOTE: This is not an official infrastructure product of the HSC collaboration.

# Stages
Available steps of the pipeline currently include:
- Acquisition of data from the HSC archive \[**get_data.py**\]
- Splitting of metadata into the different HSC fields (and optionally by photometric filter) \[**split_metadata.py**\]
- Cleaning of the raw catalogues \[**clean_catalogues.py**\]
- Applying selection cuts to define samples for clustering analysis \[**sample_selection.py**\] 
- Making maps of various quantities from the survey metadata using [decasu](https://github.com/erykoff/decasu) \[**make_maps_from_metadata.py**\]
- Create maps of the survey footprint, depth, and mask, along with maps of dust attenuation in each band and star counts, using information from the galaxy catalogues. \[**make_maps_from_catalogue.py**\]
- Create maps of the galaxy number counts and density contrast \[**make_galaxy_maps.py**\]
- Combine the maps from the individual HSC fields \[**combine_fields.py**\]
- Perform principal component analysis (PCA) on the systematics maps produced in previous stages \[**pca_systematics**\]
- Estimating redshift distributions for each galaxy sample using the 'DIR' method with COSMOS2020 ([Weaver et al. 2022](https://arxiv.org/abs/2110.13923)) photometric redshifts \[**dir_photozs.py**\]
- Compute theoretical predictions for the clustering angular power spectra using [CCL](https://github.com/LSSTDESC/CCL) \[**theory_predictions.py**\]
- Compute the clustering angular power from the data using [NaMaster](https://github.com/LSSTDESC/NaMaster) \[**compute_power_spectra.py**\]
- Plot the angular power spectra \[**plot_power_spectra.py**\]
- Convert products from this pipeline into appropriate formats for [TXPipe](https://github.com/LSSTDESC/TXPipe) \[**make_txpipe_inputs.py**\]

# Configuration
Configuration of the pipeline is done using a `yaml` file. This file will contain global settings under the heading `global`, and settings for individual stages under headings corresponding to the name of that stage. Example config files can be found in `pipelines/`. 

# Running locally
If running on a local machine, users can choose which stages to run by toggling the `True/False` settings located in **MASTER.py**.

# Running on NERSC
If running on NERSC, users must instead choose which stages to run by commenting/uncommenting the appropriate stages in **MASTER_NERSC.sh**. Executing this script will then run the appropriate job script from the **slurm_scripts** directory. 
