import os


def add_photoz(mt):
    quants = ['mean', 'mode', 'best', 'mc', 'err68_min', 'err68_max',
              'err95_min', 'err95_max']
    sto = ',\n\t'.join([f'p{mt}.photoz_{q} as pz_{q}_{mt}' for q in quants])
    return sto


def submit_job(
        fname_sql,
        output_format,
        do_preview=False,
        do_download=True,
        do_delete=True,
        output_file='query.sql',
        release='pdr3-citus-columnar'
        ):

    command = "python hscReleaseQueryDR3.py --user=tmcornish "\
        f"--release-version={release}"
    if do_preview:
        command += " -p"
    command += " -f "+output_format
    if do_delete:
        command += " -D"
    if do_download:
        if os.path.isfile(output_file):
            print("Found "+output_file)
            return
        command += " -d"
    command += " "+fname_sql
    if do_download:
        command += " > "+output_file

    print(command)
    os.system(command)


def write_frames(tablename, fname_out, dir_out='./', output_format='fits',
                 submit=False, do_download=False):

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
        submit_job(fname_out, output_format, do_download=do_download,
                   output_file=fname_job)


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
        mizuki_mstar_sfrs=False,
        apply_cuts=False,
        strict_cuts=False,
        do_download=False,
        part=None
        ):

    def add_filters(name, cat='forced', behind=False):
        if behind:
            l_str = [f'{cat}.{name}{fi}' for fi in filters]
        else:
            l_str = [f'{cat}.{fi}{name}' for fi in filters]

        sthere = ',\n\t'.join(l_str)
        return sthere

    stout_final = f'-- Run field, {fname_out}'+'\n'

    if do_download:
        fname_job = f'{dir_out}{tablename.upper()}_{fieldname.upper()}'

        if part is not None:
            fname_job += f'_part{part:02d}'
        if not strict_cuts:
            fname_job += '_shearcat'

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
        'mag_corr.corr_rmag',
        'mag_corr.corr_imag',
        'meas2.i_blendedness_abs as i_blendedness_abs',
        add_filters('a_', behind=True),
        add_filters('_psfflux_flux', cat='forced2'),
        add_filters('_psfflux_fluxerr', cat='forced2'),
        add_filters('_psfflux_flag', cat='forced2'),
        add_filters('_psfflux_mag', cat='forced2'),
        add_filters('_psfflux_magerr', cat='forced2'),
        add_filters('_cmodel_flux'),
        add_filters('_cmodel_fluxerr'),
        add_filters('_cmodel_flux'),
        add_filters('_cmodel_flag'),
        add_filters('_cmodel_mag'),
        add_filters('_cmodel_magerr'),
        add_filters('_mask_brightstar_ghost15', cat='masks'),
        add_filters('_mask_brightstar_halo', cat='masks'),
        add_filters('_mask_brightstar_blooming', cat='masks'),
        add_filters('_mask_brightstar_dip', cat='masks'),
        'masks.y_mask_brightstar_channel_stop',
    ]

    if not apply_cuts:
        stout += [
            'forced.isprimary',
            'forced.i_cmodel_flag_badcentroid',
            'meas.i_deblend_skipped',
            ]
        if strict_cuts:
            stout += [
                'forced.i_cmodel_flag_badcentroid',
                'forced2.i_sdsscentroid_flag',
                'forced.i_pixelflags_edge',
                'forced.i_pixelflags_interpolatedcenter',
                'forced.i_pixelflags_saturatedcenter',
                'forced.i_pixelflags_crcenter',
                'forced.i_pixelflags_bad',
                'forced.i_pixelflags_suspectcenter',
                'forced.i_pixelflags_clipped',
            ]
        else:
            stout += [
                add_filters('_cmodel_flag_badcentroid'),
                add_filters('_sdsscentroid_flag', cat='forced2'),
                add_filters('_pixelflags_edge'),
                add_filters('_pixelflags_interpolatedcenter'),
                add_filters('_pixelflags_saturatedcenter'),
                add_filters('_pixelflags_crcenter'),
                add_filters('_pixelflags_bad'),
                add_filters('_pixelflags_suspectcenter'),
                add_filters('_pixelflags_clipped'),
            ]

    if do_photoz:
        pzcodes = ['demp', 'dnnz', 'mizu']
        stout += [add_photoz(mt) for mt in pzcodes]
    if mizuki_mstar_sfrs:
        stout += [
            'pmizu.stellar_mass as mstar_mizu',
            'pmizu.stellar_mass_err68_min as mstar_err68_min_mizu',
            'pmizu.stellar_mass_err68_max as mstar_err68_max_mizu',
            'pmizu.sfr as sfr_mizu',
            'pmizu.sfr_err68_min as sfr_err68_min_mizu',
            'pmizu.sfr_err68_max as sfr_err68_max_mizu',
        ]
    stout = ',\n\t'.join(stout)
    stout_final += stout

    stout_final += '\nFROM\n\t'
    stout2 = [
        f'{tablename}.forced as forced',
        f'LEFT JOIN {tablename}.forced2 USING (object_id)',
        f'LEFT JOIN {tablename}.forced3 USING (object_id)',
        f'LEFT JOIN {tablename}.mag_corr USING (object_id)',
        f'LEFT JOIN {tablename}.meas meas USING (object_id)',
        f'LEFT JOIN {tablename}.meas2 meas2 USING (object_id)',
        f'LEFT JOIN {tablename}.masks masks USING (object_id)'
    ]
    if do_photoz:
        pzcodes_full = ['demp', 'dnnz', 'mizuki']
        stout2 += [f'LEFT JOIN {tablename}.photoz_{mt_full} p{mt} '
                   'USING (object_id)'
                   for mt, mt_full in zip(pzcodes, pzcodes_full)]
    elif mizuki_mstar_sfrs:
        stout2 += [f'LEFT JOIN {tablename}.photoz_mizuki pmizu '
                   'USING (object_id)'
                   for mt, mt_full in zip(pzcodes, pzcodes_full)]
    stout2 = '\n\t'.join(stout2)
    stout_final += stout2

    if apply_cuts or do_field or (ra_range is not None) or \
            (dec_range is not None):
        stout_final += '\nWHERE\n\t'
        stout3 = []
        if apply_cuts:
            stout3 += [
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

        stout_final += stout3

    stout_final += '\n;\n'

    with open(fname_out, 'w') as f:
        f.write(stout_final)

    if submit:
        submit_job(fname_out, output_format, do_download=do_download,
                   output_file=fname_job)
