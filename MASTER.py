#####################################################################################################
# Master script for running different stages of the HSC-SSP DR3 analysis.
# If running on glamdring, run 'source /mnt/zfsusers/tcornish/venvs/lsst/bin/activate' before this.
#
# TODO: incorporate option to run via SLURM somehow
#####################################################################################################

# Import necessary packages/modules
import os, sys
import output_utils as opu


##################
#### SETTINGS ####
##################

#toggle `switches' for determining which scripts to run
get_data = False		#run data acquisition script
split_meta = False		#splits metadata by global field
clean_cats = False		#apply various cuts to clean the catalogues
catbased_maps = False	#make maps for various quantities using the catalogue
metadata_maps = False	#make maps for various quantities using the frame metadata (uses decasu)
pca_systs = True				#perform PCA to potentially reduce the number of maps being deprojected
galaxy_maps = False		#make galaxy count and density maps in tomographic bins
power_spectra = True	#compute power spectra
txpipe_inputs = False	#collects all relevant files and converts them into TXPipe-compatible formats


####################

settings = [
	get_data,
	split_meta,
	clean_cats,
	catbased_maps,
	metadata_maps,
	pca_systs,
	galaxy_maps,
	power_spectra,
	txpipe_inputs
	]

proc = [
	'Downloading data from HSC database',
	'Splitting metadata by field',
	'Cleaning catalogues',
	'Making maps from catalogue data',
	'Making maps from frame metadata',
	'Performing PCA',
	'Making galaxy count and density maps in z bins',
	'Computing power spectra',
	'Making TXPipe-compatible inputs'
	]

run_str = [
	'cd data_query/ && python -u get_data.py; cd ..',
	'python -u split_metadata.py',
	'python -u clean_catalogues.py',
	'python -u make_maps_from_catalogue.py',
	'python -u make_maps_from_metadata.py',
	'python -u pca_systematics.py',
	'python -u make_galaxy_maps.py',
	'python -u compute_power_spectra.py',
	'python -u make_txpipe_inputs.py'
	]



print(opu.colour_string(opu.string_important('PROCESSES TO RUN')+'\n', 'cyan'))
setting_str = []
for se, pn in zip(settings, proc):
	if se:
		setting_str.append(pn)
print('\n'.join(setting_str)+'\n')



#########################
#### RUNNING SCRIPTS ####
#########################

for se, pn, rs in zip(settings, proc, run_str):
	if se:
		print(opu.colour_string(opu.string_important(pn.upper())+'\n', 'orange')+'\n')
		os.system(rs)