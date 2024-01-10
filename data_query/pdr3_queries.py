import numpy as np
import os
import sys
import contextlib


def add_photoz(mt):
	quants = ['mean', 'mode', 'best', 'mc']
	sto = ',\n\t'.join([f'p{mt}.photoz_{q} as pz_{q}_{mt}' for q in quants])
	return sto


def submit_job(fname_sql,output_format,do_preview=False,do_download=True,do_delete=True,output_file='query.sql',release='pdr3') :

	command=f"python hscReleaseQueryDR3.py --user=tmcornish --release-version=pdr3-citus-columnar"
	if do_preview :
		command+=" -p"
	command+=" -f "+output_format
	if do_delete:
		command+=" -D"
	if do_download :
		if os.path.isfile(output_file) :
			print("Found "+output_file)
			return
		command+=" -d"
	command+=" "+fname_sql
	if do_download :
		command+=" > "+output_file

	print(command)
	os.system(command)


def write_frames(
	tablename,
	fname_out,
	dir_out='./',
	output_format='fits',
	submit=False,
	do_download=False
	):
	
	stout = [
		'SELECT *',
		f'FROM {tablename}.frame',
		';']
	stout = '\n'.join(stout)

	if do_download:
		fname_job = f'{dir_out}{tablename.upper()}_frames.{output_format}'
	else:
		fname_job = None

	with open(fname_out, 'w') as f:
		f.write(stout)

	if submit:
		submit_job(fname_out, output_format, do_download=do_download, output_file=fname_job)



def write_fieldsearch(
	tablename,
	fieldname,
	fname_out,
	dir_out='./',
	output_format='fits',
	filters=['g', 'r', 'i', 'z', 'y'],
	do_field=True,
	submit=False,
	ra_range=None,
	dec_range=None,
	do_photoz=False,
	strict_cuts=False,
	do_download=False,
	part=None
	):
	
	def add_filters(name, cat='forced', behind=False):
		if behind:
			l = [f'{cat}.{name}{fi}' for fi in filters]
		else:
			l = [f'{cat}.{fi}{name}' for fi in filters]

		sthere = ',\n\t'.join(l)
		return sthere

	stout_final = f'-- Run field, {fname_out}'+'\n'
	#stout = [f'-- Run field, {fname_out}'+'\n']

	if do_download:
		fname_job = f'{dir_out}{tablename.upper()}_{fieldname.upper()}'

		if part is not None:
			fname_job += f'_part{part}'
		if not strict_cuts:
			fname_job += f'_shearcat'

		fname_job += f'_forced.{output_format}'

	stout = [
		'SELECT object_id',
		'forced.ra as ra',
		'forced.dec as dec',
		'forced.tract as tract',
		'forced.patch as patch',
		add_filters('merge_peak_', behind=True),
		add_filters('_inputcount_value'),
		'forced.i_pixelflags_bright_objectcenter',
		'forced.i_pixelflags_bright_object',
		'forced.i_extendedness_value',
		'meas2.i_blendedness_abs as i_blendedness_abs',
		add_filters('a_', behind=True),
		add_filters('_psfflux_flux', cat='forced2'),
		add_filters('_psfflux_fluxerr', cat='forced2'),
		add_filters('_psfflux_flag', cat='forced2'),
		add_filters('_psfflux_mag', cat='forced2'),
		add_filters('_psfflux_magerr', cat='forced2'),
		add_filters('_apertureflux_10_flux', cat='forced3'),
		add_filters('_apertureflux_10_fluxerr', cat='forced3'),
		add_filters('_apertureflux_10_flag', cat='forced3'),
		add_filters('_apertureflux_10_mag', cat='forced3'),
		add_filters('_apertureflux_10_magerr', cat='forced3'),
		add_filters('_cmodel_flux'),
		add_filters('_cmodel_fluxerr'),
		add_filters('_cmodel_flux'),
		add_filters('_cmodel_flag'),
		add_filters('_cmodel_mag'),
		add_filters('_cmodel_magerr'),
		'masks.i_mask_brightstar_ghost15',
		'masks.i_mask_brightstar_halo',
		'masks.i_mask_brightstar_blooming',
	]
	if do_photoz:
		pzcodes = ['demp', 'dnnz', 'mizu']
		stout += [add_photoz(mt) for mt in pzcodes]
	stout = ',\n\t'.join(stout)
	stout_final += stout

	stout_final += '\nFROM\n\t'
	stout2 = [
		f'{tablename}.forced as forced',
		f'LEFT JOIN {tablename}.forced2 USING (object_id)',
		f'LEFT JOIN {tablename}.forced3 USING (object_id)',
		f'LEFT JOIN {tablename}.meas meas USING (object_id)',
		f'LEFT JOIN {tablename}.meas2 meas2 USING (object_id)',
		f'LEFT JOIN {tablename}.masks masks USING (object_id)'
	]
	if do_photoz:
		pzcodes_full = ['demp', 'dnnz', 'mizuki']
		stout2 += [f'LEFT JOIN {tablename}.photoz_{mt_full} p{mt} USING (object_id)' for mt,mt_full in zip(pzcodes,pzcodes_full)]
	stout2 = '\n\t'.join(stout2)
	stout_final += stout2

	stout_final += '\nWHERE\n\t'
	stout3 = [
		'forced.isprimary=True',
		'forced.i_cmodel_flag_badcentroid=False',
		'forced2.i_sdsscentroid_flag=False',
		'forced.i_pixelflags_edge=False',
		'forced.i_pixelflags_interpolatedcenter=False',
		'forced.i_pixelflags_saturatedcenter=False',
		'forced.i_pixelflags_crcenter=False',
		'forced.i_pixelflags_bad=False',
		'forced.i_pixelflags_suspectcenter=False',
		'forced.i_pixelflags_clipped=False',
		'meas.i_deblend_skipped=False'
		]
	if strict_cuts:
		stout3 += [
		'forced2.g_sdsscentroid_flag=False',
		'forced2.r_sdsscentroid_flag=False',
		'forced2.z_sdsscentroid_flag=False',
		'forced2.y_sdsscentroid_flag=False',
		'forced.g_cmodel_flag=False',
		'forced.r_cmodel_flag=False',
		'forced.i_cmodel_flag=False',
		'forced.z_cmodel_flag=False',
		'forced.y_cmodel_flag=False',
		'forced2.g_psfflux_flag=False',
		'forced2.r_psfflux_flag=False',
		'forced2.i_psfflux_flag=False',
		'forced2.z_psfflux_flag=False',
		'forced2.y_psfflux_flag=False',
		'forced.g_pixelflags_edge=False',
		'forced.r_pixelflags_edge=False',
		'forced.z_pixelflags_edge=False',
		'forced.y_pixelflags_edge=False',
		'forced.g_pixelflags_interpolatedcenter=False',
		'forced.r_pixelflags_interpolatedcenter=False',
		'forced.z_pixelflags_interpolatedcenter=False',
		'forced.y_pixelflags_interpolatedcenter=False',
		'forced.g_pixelflags_saturatedcenter=False',
		'forced.r_pixelflags_saturatedcenter=False',
		'forced.z_pixelflags_saturatedcenter=False',
		'forced.y_pixelflags_saturatedcenter=False',
		'forced.g_pixelflags_crcenter=False',
		'forced.r_pixelflags_crcenter=False',
		'forced.z_pixelflags_crcenter=False',
		'forced.y_pixelflags_crcenter=False',
		'forced.g_pixelflags_bad=False',
		'forced.r_pixelflags_bad=False',
		'forced.z_pixelflags_bad=False',
		'forced.y_pixelflags_bad=False'
		]
	if do_field:
		stout3 += [f"forced.field='{fieldname}'"]
	
	if ra_range is not None:
		stout3 += [
			f'forced.ra>={ra_range[0]:.3f}',
			f'forced.ra<{ra_range[1]:.3f}'
			]
	if dec_range is not None:
		stout3 += [
			f'forced.dec>={dec_range[0]:.3f}',
			f'forced.dec<{dec_range[1]:.3f}'
			]
	
	stout3 = ' and\n\t'.join(stout3)

	#stout_final = '\n'.join([stout, stout2, stout3]) + '\n;\n'
	stout_final += stout3 + '\n;\n'

	with open(fname_out, 'w') as f:
		f.write(stout_final)

	if submit:
		submit_job(fname_out, output_format, do_download=do_download, output_file=fname_job)

#write_fieldsearch('pdr3_wide', 'hectomap', 'HECTOMAP_query_test.sql', dir_out='/Users/thomascornish/Desktop/LSST_clustering/Data/HSC_DR3/', submit=False, strict_cuts=True, do_photoz=True)


