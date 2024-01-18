#####################################################################################################
# Master script for running different stages of the HSC-SSP DR3 analysis.
# If running on glamdring, run 'source /mnt/zfsusers/tcornish/venvs/lsst/bin/activate' before this.
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
clean_cats = True		#apply various cuts to clean the catalogues
split_by_pixel = False	#split the catalogues by pixel 
catbased_maps = True	#make maps for various quantities using the catalogue
metadata_maps = False	#make maps for various quantities using the frame metadata (uses decasu)
galaxy_maps = True		#make galaxy count and density maps in tomographic bins


####################

settings = [
	get_data,
	split_meta,
	clean_cats,
	split_by_pixel,
	catbased_maps,
	metadata_maps,
	galaxy_maps
	]

proc = [
	'Downloading data from HSC database',
	'Splitting metadata by field',
	'Cleaning catalogues',
	'Splitting data by pixel',
	'Making maps from catalogue data',
	'Making maps from frame metadata',
	'Making galaxy count and density maps in z bins'
	]

run_str = [
	'cd data_query/ && python get_data.py; cd ..',
	'python split_metadata.py',
	'python clean_catalogues.py',
	'python split_data_by_pixel.py',
	'python make_maps_from_catalogue.py',
	'python make_maps_from_metadata.py',
	'python make_galaxy_maps.py'
	]

#TODO: change to output bash script which runs relevant scripts
#TODO: add option for above bash script to be formatted for SLURM (i.e. glamdring)

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