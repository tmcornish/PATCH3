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
clean_cats = False		#apply various cuts to clean the catalogues
catbased_maps = False	#make maps for various quantities using the catalogue
metadata_maps = True	#make maps for various quantities using the frame metadata (uses decasu)
galaxy_maps = False		#make galaxy count and density maps in tomographic bins


####################

settings = [
	get_data,
	clean_cats,
	catbased_maps,
	metadata_maps,
	galaxy_maps
	]

proc = [
	'Downloading data from HSC database',
	'Cleaning catalogues',
	'Making maps from catalogue data',
	'Making maps from frame metadata',
	'Making galaxy count and density maps in z bins'
	]

run_str = [
	'cd data_query/ && python get_data.py; cd ..',
	'python clean_catalogues.py',
	'python make_maps_from_catalogue.py',
	'python make_maps_from_metadata.py',
	'python make_galaxy_maps.py'
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