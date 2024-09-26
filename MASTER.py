#####################################################################################################
# Master script for running different stages of the HSC-SSP DR3 analysis.
#####################################################################################################

# Import necessary packages/modules
import os, sys
import output_utils as opu


##################
#### SETTINGS ####
##################

#toggle `switches' for determining which scripts to run
get_data = False		#run data acquisition script
split_meta = False		#splits metadata by field
clean_cats = False		#apply various cuts to clean the catalogues
metadata_maps = False	#make maps for various quantities using the frame metadata (uses decasu)
catbased_maps = False	#make maps for various quantities using the catalogue
galaxy_maps = False		#make galaxy count and density maps in tomographic bins
combine_fields = False	#combine maps from all fields
pca_systs = False		#perform PCA to potentially reduce the number of maps being deprojected
theory_cells = False		#compute theoretical angular power spectra
power_spectra = False	#compute power spectra
make_sacc = False		#consolidate c_ell info into a SACC file
plot_cells = True		#plot the power spectra
txpipe_inputs = False	#collects all relevant files and converts them into TXPipe-compatible formats


####################

settings = [
	get_data,
	split_meta,
	clean_cats,
	metadata_maps,
	catbased_maps,
	galaxy_maps,
	combine_fields,
	pca_systs,
	theory_cells,
	power_spectra,
	make_sacc,
	plot_cells,
	txpipe_inputs
	]

proc = [
	'Downloading data from HSC database',
	'Splitting metadata by field',
	'Cleaning catalogues',
	'Making maps from frame metadata',
	'Making maps from catalogue data',
	'Making galaxy count and density maps in z bins',
	'Combining maps from all fields',
	'Performing PCA',
	'Computing theoretical power spectra',
	'Computing power spectra',
	'Creating SACC file',
	'Plotting power spectra',
	'Making TXPipe-compatible inputs'
	]

run_str = [
	'cd data_query/ && python -u get_data.py; cd ..',
	'python -u split_metadata.py',
	'python -u clean_catalogues.py',
	'python -u make_maps_from_metadata.py',
	'python -u make_maps_from_catalogue.py',
	'python -u make_galaxy_maps.py',
	'python -u combine_fields.py',
	'python -u pca_systematics.py',
	'python -u theory_predictions.py',
	'python -u compute_power_spectra.py',
	'python -u make_sacc_files.py',
	'python -u plot_power_spectra.py',
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