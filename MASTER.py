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
split_meta = True		#splits metadata by global field
clean_cats = True		#apply various cuts to clean the catalogues
#split_by_pixel = False	#split the catalogues by pixel 
catbased_maps = True	#make maps for various quantities using the catalogue
#metadata_maps = False	#make maps for various quantities using the frame metadata (uses decasu)
galaxy_maps = True		#make galaxy count and density maps in tomographic bins
power_spectra = False	#compute power spectra


####################

settings = [
	get_data,
	split_meta,
	clean_cats,
	#split_by_pixel,
	catbased_maps,
	#metadata_maps,
	galaxy_maps,
	power_spectra
	]

proc = [
	'Downloading data from HSC database',
	'Splitting metadata by field',
	'Cleaning catalogues',
	#'Splitting data by pixel',
	'Making maps from catalogue data',
	#'Making maps from frame metadata',
	'Making galaxy count and density maps in z bins',
	'Computing power spectra'
	]

run_str = [
	'cd data_query/ && python -u get_data.py; cd ..',
	'python -u split_metadata.py',
	'python -u clean_catalogues.py',
	#'python -u split_data_by_pixel.py',
	'python -u make_maps_from_catalogue.py',
	#'python -u make_maps_from_metadata.py',
	'python -u make_galaxy_maps.py'
	'python -u compute_power_spectra.py'
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