##########################################################################
# - Applies various cuts to each of the catalogues downloaded from the HSC
#   database.
# - Combines all cleaned catalogues belonging to the same field.
##########################################################################

import os
import sys
from astropy.table import Table, vstack
import numpy as np
import pandas as pd
import glob
import h5py
from configuration import PipelineConfig as PC
from output_utils import colour_string, write_output_hdf, \
    h5py_dataset_iterator, error_message

# SETTINGS #
config_file = sys.argv[1]
cf = PC(config_file, stage='cleanCats')


###################
#    FUNCTIONS    #
###################

def basic_clean(t):
    '''
    Performs a basic clean of the input catalogue using pre-existing flags.

    Parameters
    ----------
    t: astropy.table.Table
        Input catalogue.

    Returns
    -------
    t: astropy.table.Table
        Catalogue after applying the cuts.
    '''

    sel = np.ones(len(t), dtype=bool)

    # Create an empty list to which column names with the 'is_null' suffix
    # will be appended
    isnull_names = []

    # See if file exists containing a list of required columns
    if os.path.exists(cf.auxfiles.required_cols):
        with open(cf.auxfiles.required_cols) as file:
            req_cols = file.readlines()
        # Remove newline escape sequences and delete any commented lines
        req_cols = [c.replace('\n', '') for c in req_cols
                    if not c.startswith('#')]
    else:
        req_cols = t.colnames
    # Cycle through column names
    for key in t.colnames:
        if key.__contains__('isnull') and key in req_cols:
            sel[t[key]] = 0
            isnull_names.append(key)
        else:
            # Remove NaNs UNLESS they are photo-zs, SFRs, stellar masses
            if not key.startswith('pz_') and not key.startswith('corr_') \
                    and not key.startswith('mstar_') \
                    and not key.startswith('sfr_'):
                sel[np.isnan(t[key])] = 0

    # Also remove any sources that are not primary detections
    sel *= t['isprimary']

    # If told to apply 'main' and 'strict'  cuts at this stage, do so
    if len(cf.remove_if_flagged) > 0:
        import itertools
        to_remove = cf.combine_flags(
                        t,
                        list(itertools.chain(*[cf.flags[k] for k in cf.flags
                                               if k in cf.remove_if_flagged])),
                        combine_type='or'
        )
        sel *= ~to_remove

    t.remove_columns(isnull_names)
    t = t[sel]

    # Correct r/i-band magnitudes to r2/i2 using the appropriate corrections
    if cf.correct_ri:
        for b in ['r', 'i']:
            corr_mask = ~np.isnan(t[f'corr_{b}mag'])
            t[f'{b}_cmodel_mag'][corr_mask] += t[f'corr_{b}mag'][corr_mask]

    return t


def photom_cuts(t):
    '''
    Performs a photometric cuts on the catalogue.

    Parameters
    ----------
    t: astropy.table.Table
        Input catalogue.

    Returns
    -------
    t: astropy.table.Table
        Catalogue after applying the cuts.
    '''

    # Magnitude cut in primary band
    sel_maglim = np.ones(len(t), dtype=bool)
    maglim_mask = (t[f'{cf.bands.primary}_cmodel_mag']
                   - t[f'a_{cf.bands.primary}']) > cf.depth_cut
    sel_maglim[maglim_mask] = False

    # Blending cut
    sel_blend = np.ones(len(t), dtype=bool)
    blend_mask = t[
        f'{cf.bands.primary}_blendedness_abs'
        ] >= 10 ** cf.log_blend_cut
    sel_blend[blend_mask] = False

    # S/N cut in primary band
    sel_sn_pri = np.ones(len(t), dtype=bool)
    sn_pri_mask = (t[f'{cf.bands.primary}_cmodel_flux']
                   / t[f'{cf.bands.primary}_cmodel_fluxerr']) < cf.sn_pri
    sel_sn_pri[sn_pri_mask] = False

    # S/N cut in all other bands
    sel_sn_sec = []
    for b in cf.bands.secondary:
        sel_sn_now = np.ones(len(t), dtype=bool)
        sn_sec_mask = (t[f'{b}_cmodel_flux']
                       / t[f'{b}_cmodel_fluxerr']) < cf.sn_sec
        sel_sn_now[sn_sec_mask] = False
        sel_sn_sec.append(sel_sn_now)
    # Ensure that source is above threshold in at least two of these bands
    sel_sn_sec = np.sum(sel_sn_sec, axis=0) >= 2

    # Apply the cuts
    t = t[sel_maglim * sel_blend * sel_sn_pri * sel_sn_sec]

    return t


def gal_cut(t):
    '''
    Flags sources as stars or galaxies based on an 'extendedness' parameter.

    Parameters
    ----------
    t: astropy.table.Table
        Input catalogue.

    Returns
    -------
    t_gals: astropy.table.Table
        Catalogue containing only the sources identified as galaxies.

    t_stars: astropy.table.Table
        Catalogue containing only the sources identified as stars.
    '''

    # Identify stars via their 'extendedness' in the primary band
    star_mask = t[f'{cf.bands.primary}_extendedness_value'] == 0.

    # Split the catalogue into stars and galaxies
    t_stars = t[star_mask]
    t_gals = t[~star_mask]

    return t_gals, t_stars


def flag_stars(t):
    '''
    Adds a column to the provided table flagging all stars.

    Parameters
    ----------
    t: astropy.table.Table
        Input catalogue.

    '''

    # Identify stars via their 'extendedness' in the primary band
    star_mask = t[f'{cf.bands.primary}_extendedness_value'] == 0.
    # Add this as a column to the Table
    t['is_star'] = star_mask


def get_data_suffix():
    '''
    Determines the suffix used when saving the downloaded raw data.
    '''
    suffix = ''
    if not cf.strict_cuts:
        suffix += '_shearcat'
    suffix += '_forced'
    return suffix


#######################################################
#                  START OF SCRIPT                    #
#######################################################

# Descriptors for each sample during cleaning
catalogues = [
    'raw',
    'basic_cleaned',
    'fully_cleaned',
    'galaxies',
    'stars'
]

# Get a dictionary of all fields being analysed and their respective subfields
fields = cf.fields

# Cycle through each global field
for g in fields:
    cf.fields = [g]
    PATH_G = f'{cf.paths.out}{g}'
    # Counters for keeping track of the cleaning stages in the whole field
    l_init_fd = 0 		# Counter for number of sources in raw data
    l_bc_fd = 0			# Counter for number of sources in basic-cleaned data
    l_final_fd = 0 		# Counter for number of sources in final catalogue
    l_gals_fd = 0		# Counter for number of galaxies in final catalogue
    l_stars_fd = 0		# Counter for number of stars in final catalogue
    print(colour_string(g.upper(), 'orange'))
    # Cycle through each of the subfields
    for fd in cf.get_subfields():
        print(colour_string(fd, 'purple'))

        # Create output directory for this field
        OUT = f'{PATH_G}/{fd}'
        print(f'Output directory: {OUT}')
        if not os.path.exists(OUT):
            os.system(f'mkdir -p {OUT}')

        # Filenames to be given to the HDF format output files
        hdf_basic = f'{OUT}/{cf.cats.basic}'
        hdf_full = f'{OUT}/{cf.cats.main}'
        # See if the field has been split into multiple parts
        fname = f'{cf.paths.data}{cf.dr.upper()}_{fd.upper()}'\
                f'{get_data_suffix()}.fits'
        # Initially enable 'write' mode for output files
        mode = 'w'
        if os.path.exists(fname):
            data_all = Table.read(fname, format='fits')
            l_init = len(data_all)
            # Apply basic clean and write to HDF file
            print('Applying basic clean...')
            data_all = basic_clean(data_all)
            flag_stars(data_all)
            l_bc = len(data_all)
            write_output_hdf(data_all, hdf_basic, mode=mode,
                             group='photometry')
            # Apply photometric cuts and write to HDF file
            print('Applying photometric cuts...')
            data_all = photom_cuts(data_all)
            write_output_hdf(data_all, hdf_full, mode=mode, group='photometry')
            l_final = len(data_all)
        else:
            # See if catalogues exist for separate parts of the field
            parts = sorted(glob.glob(f'{cf.paths.data}{cf.dr.upper()}'
                                     f'_{fd.upper()}_part??{get_data_suffix()}'
                                     '.fits'))
            if len(parts) >= 1:
                # Set up a list to contain data from all catalogues associated
                # with this field
                data_all = []
                l_init = 0		# Counter for number of sources in raw data
                l_bc = 0		# Counter for number of sources in basic-cleaned data
                l_final = 0		# Counter for number of sources in final catalogue
                # Cycle through each catalogue
                for i, cat in enumerate(parts):
                    print(f'Cleaning part {i+1}...')
                    data = Table.read(cat, format='fits')
                    l_init += len(data)
                    # Apply basic clean and write to HDF file
                    print('Applying basic clean...')
                    data = basic_clean(data)
                    flag_stars(data)
                    l_bc += len(data)
                    write_output_hdf(data, hdf_basic, mode=mode,
                                     group='photometry')
                    # Apply photometric cuts and write to HDF file
                    print('Applying photometric cuts...')
                    data = photom_cuts(data)
                    write_output_hdf(data, hdf_full, mode=mode,
                                     group='photometry')
                    data_all.append(data)
                    l_final += len(data)
                    try:
                        assert mode == 'a'
                    except AssertionError:
                        mode = 'a'
                # Stack the data from each part
                data_all = vstack(data_all)

            else:
                error_message(cf.__name__,
                              f'No catalogues found for field {fd.upper()}.')
                continue

        print(colour_string(f'Began with {l_init} sources.', 'green'))
        print(colour_string(f'{l_bc} remained after basic cleaning.', 'green'))
        print(
            colour_string(f'{l_final} sources remaining after full cleaning.',
                          'green')
            )

        # Split catalogue into stars and galaxies
        data_gals, data_stars = gal_cut(data_all)
        l_gals = len(data_gals)
        l_stars = len(data_stars)
        print(colour_string(f'{l_gals} galaxies; {l_stars} stars.', 'green'))

        # Write the catalogues to output files
        print('Writing outputs...')
        hdf_stars = f'{OUT}/{cf.cats.stars}'
        write_output_hdf(data_gals, hdf_full, mode='w', group='photometry')
        write_output_hdf(data_stars, hdf_stars, mode='w', group='photometry')

        # Add to the relevant counters
        l_init_fd += l_init
        l_bc_fd += l_bc
        l_final_fd += l_final
        l_gals_fd += l_gals
        l_stars_fd += l_stars

    print(colour_string(f'SUMMARY: {g.upper()}', 'orange'))
    print(colour_string(f'Began with {l_init_fd} sources.', 'green'))
    print(colour_string(f'{l_bc_fd} remained after basic cleaning.', 'green'))
    print(colour_string(f'{l_final_fd} sources remaining after full cleaning.',
                        'green'))
    print(colour_string(f'{l_gals_fd} galaxies; {l_stars_fd} stars.', 'green'))

    # Write summary of each stage to file
    summary_data = [
        l_init_fd,
        l_bc_fd,
        l_final_fd,
        l_gals_fd,
        l_stars_fd
        ]
    summary = {'Catalogue': catalogues, 'N_sources': summary_data}
    summary = pd.DataFrame(data=summary)
    summary.to_csv(f'{PATH_G}/{cf.clean_summary_file}', sep='\t', index=False)

print('Consolidating catalogues from subfields...')
cats = [cf.cats.basic, cf.cats.main, cf.cats.stars]
for g in fields:
    cf.fields = [g]
    print(colour_string(g.upper(), 'orange'))
    # Cycle through the catalogue types
    for cat in cats:
        print(colour_string(cat, 'cyan'))
        fname = f'{cf.paths.out}{g}/{cat}'
        with h5py.File(fname, 'w') as fmain:
            for fd in cf.get_subfields():
                print(f'Adding data from subfield {fd}...')
                cat_now = f'{cf.paths.out}{g}/{fd}/{cat}'
                with h5py.File(cat_now, 'r') as fnow:
                    # Check the current structure of the main file
                    paths_current = [p for p, _
                                     in h5py_dataset_iterator(fmain)]
                    # Iterate through each branch of the file tree
                    for (path, dset) in h5py_dataset_iterator(fnow):
                        N = len(dset)
                        dt = dset.dtype
                        if path in paths_current:
                            dset_main = fmain[path]
                            dset_main.resize((len(dset_main)+N,))
                            dset_main[-N:] = dset[:]
                        else:
                            dset = fmain.create_dataset(path, shape=(N,),
                                                        data=dset[:],
                                                        maxshape=(None,),
                                                        dtype=dt)
                if cf.remove_intermediate:
                    os.system(f'rm -f {cat_now}')
