###########################################################################
# Master script for running different stages of the HSC-SSP DR3 analysis.
###########################################################################

# Import necessary packages/modules
import os
import sys
import output_utils as opu

# Retrieve the path of the pipeline config file
config_file = sys.argv[1]

##################
#    SETTINGS    #
##################

# Toggle `switches' for determining which scripts to run
get_data = False        # Run data acquisition script
split_meta = False      # Splits metadata by field
clean_cats = False      # Apply various cuts to clean the catalogues
select_samples = True   # Apply sample selection in clean catalogues
metadata_maps = False   # Make maps from frame metadata (uses decasu)
catbased_maps = False   # Make maps from quantities in the catalogues
galaxy_maps = False     # Make galaxy count and density maps for each sample
combine_fields = False  # Combine maps from all fields
pca_systs = False       # Perform PCA on systematics templates
dir_photozs = False     # Use DIR to compute n(z) distributions
theory_cells = False    # Compute theoretical angular power spectra
power_spectra = False   # Compute power spectra
fit_hods = False        # Fit HODs to the computed angular power spectra
plot_cells = False      # Plot the power spectra
txpipe_inputs = False   # Convert outputs into TXPipe-compatible formats


####################

settings = [
    get_data,
    split_meta,
    clean_cats,
    select_samples,
    metadata_maps,
    catbased_maps,
    galaxy_maps,
    combine_fields,
    pca_systs,
    dir_photozs,
    theory_cells,
    power_spectra,
    fit_hods,
    plot_cells,
    txpipe_inputs
    ]

proc = [
    'Downloading data from HSC database',
    'Splitting metadata by field',
    'Cleaning catalogues',
    'Applying sample selection',
    'Making maps from frame metadata',
    'Making maps from catalogue data',
    'Making galaxy count and density maps in z bins',
    'Combining maps from all fields',
    'Performing PCA',
    'Computing n(z) distributions using DIR',
    'Computing theoretical power spectra',
    'Computing power spectra',
    'Fitting HOD model',
    'Plotting power spectra',
    'Making TXPipe-compatible inputs'
    ]

run_str = [
    f'python -u -m data_query.get_data {config_file}',
    f'python -u split_metadata.py {config_file}',
    f'python -u clean_catalogues.py {config_file}',
    f'python -u sample_selection.py {config_file}',
    f'python -u make_maps_from_metadata.py {config_file}',
    f'python -u make_maps_from_catalogue.py {config_file}',
    f'python -u make_galaxy_maps.py {config_file}',
    f'python -u combine_fields.py {config_file}',
    f'python -u pca_systematics.py {config_file}',
    f'python -u dir_photozs.py {config_file}',
    f'python -u theory_predictions.py {config_file}',
    f'python -u compute_power_spectra.py {config_file}',
    f'python -u fit_hods.py {config_file}',
    f'python -u plot_power_spectra.py {config_file}',
    f'python -u make_txpipe_inputs.py {config_file}'
    ]

print(opu.colour_string(opu.string_important('PROCESSES TO RUN')+'\n', 'cyan'))
setting_str = []
for se, pn in zip(settings, proc):
    if se:
        setting_str.append(pn)
print('\n'.join(setting_str)+'\n')

#########################
#    RUNNING SCRIPTS    #
#########################

for se, pn, rs in zip(settings, proc, run_str):
    if se:
        print(opu.colour_string(
              opu.string_important(pn.upper())+'\n', 'orange')
              + '\n')
        os.system(rs)
