\*\*\* **WORK IN PROGRESS** \*\*\*
Pipeline for performing Fourier-space galaxy clustering analyses using Data Release 3 (DR3) from the Hyper Suprime-Cam Subaru Strategic Program (HSC-SSP).

# Stages
Available steps of the pipeline currently include:
- Acquisition of data from the HSC archive \[**get_data.py**\]
- Splitting of metadata into the different HSC fields (and optionally by photometric filter) \[**split_metadata.py**\]
- Cleaning of the raw catalogues \[**clean_catalogues.py**\]
- Making maps of various quantities from the survey metadata using [decasu](https://github.com/erykoff/decasu) \[**make_maps_from_metadata.py**\]
- Create maps of the survey footprint, depth, and mask, along with maps of dust attenuation in each band and star counts, using information from the galaxy catalogues. \[**make_maps_from_catalogue.py**\]
- Create maps of the galaxy number counts and density contrast \[**make_galaxy_maps.py**\]
- Combine the maps from the individual HSC fields \[**combine_fields.py**\]
- Perform principal component analysis (PCA) on the systematics maps produced in previous stages \[**pca_systematics**\]
- Compute theoretical predictions for the clustering angular power spectra using [CCL](https://github.com/LSSTDESC/CCL) \[**theory_predictions.py**\]
- Compute the clustering angular power from the data using [NaMaster](https://github.com/LSSTDESC/NaMaster) \[**compute_power_spectra.py**\]
- Consolidate the angular power spectra information into [Sacc](https://github.com/LSSTDESC/sacc) files \[**make_sacc_files.py**\]
- Plot the angular power spectra \[**plot_power_spectra.py**\]
- Convert products from this pipeline into appropriate formats for [TXPipe](https://github.com/LSSTDESC/TXPipe)

# Configuration
Each stage of the pipeline has a corresponding `class` in **config.py**, wherein one can find and modify the settings for that stage. Global settings for the pipeline can be found in the `cf_global` class, which is located near the top of the file. All other configuration classes inherit settings from `cf_global`. 

# Running locally
If running on a local machine, users can choose which stages to run by toggling the `True/False` settings located in **MASTER.py**.

# Running on NERSC
If running on NERSC, users must instead choose which stages to run by commenting/uncommenting the appropriate stages in **MASTER_NERSC.sh**. Executing this script will then run the appropriate job script from the **slurm_scripts** directory. 
