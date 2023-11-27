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
clean_cats = True		#apply various cuts to clean the catalogues
make_maps = False		#make maps for various quantities


####################

settings = [
	get_data,
	clean_cats,
	make_maps
	]

proc = [
	'Downloading data from HSC database',
	'Cleaning catalogues',
	'Making maps'
	]

run_str = [
	'cd data_query/ && python get_data.py; cd ..',
	'python clean_catalogues.py',
	'python make_maps.py'
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