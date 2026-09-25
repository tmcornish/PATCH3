##############################################################################
# Classes for querying the HSC database via SQL.
##############################################################################

import os
from stage import baseStage


############################
#          STAGES          #
############################

class queryBase(baseStage):
    '''
    Base class for generating queries; should not be used directly.
    '''
    def __init__(
        self,
        config_file,
    ):
        super().__init__(config_file)
        self.queries = []
        self.outfiles = []
        self.stout = {}
        self.stout_str = ''
        # Output directory fro SQL query files (NOT the data itself)
        self.path_queries = self.config.paths.out + 'sql_queries/'
        # Base directory for the data
        self.path_out = self.config.paths.data
        # Mapping from different flux types to the tables to which they belong
        aper_sizes = [10, 15, 20, 30, 40, 57, 84, 118, 168, 235]
        seeings = [0, 1, 2, 3]
        conv_sizes = [11, 15, 20, 'kron']
        self.flux_tables = {
            'cmodel': 'forced',
            'gaussianflux': 'forced2',
            'psfflux': 'forced2',
            'kronflux': 'forced2',
            'sdssshape': 'forced2',
            'undeblended_psfflux': 'forced2',
            'undeblended_kronflux': 'forced2'
        } | {
            f'apertureflux_{i}': 'forced3' for i in aper_sizes
        } | {
            f'undeblended_apertureflux_{i}': 'forced3' for i in aper_sizes
        } | {
            f'convolvedflux_{s}_{c}': 'forced4'
            for s in seeings for c in conv_sizes
        } | {
            f'undeblended_convolvedflux_{s}_{c}': 'forced5'
            for s in seeings for c in conv_sizes
        }

    def _assemble_query(self):
        '''
        Assembles all components of the query into a single string.
        '''
        stout_cols = 'SELECT ' + ',\n\t'.join(self.stout['cols'])
        stout_from = 'FROM ' + '\n\tLEFT JOIN '.join(self.stout['from'])
        stout_conds = 'WHERE ' + ' AND \n\t'.join(self.stout['conds'])
        self.stout_str = '\n'.join(
            [
                stout_cols,
                stout_from,
                stout_conds
            ]
        ) + '\n;'

    def write_sql(self):
        '''
        Generates a query and writes it to an SQL file.

        Here this is just a placeholder; it is properly defined for
        each subclass.
        '''
        print('WARNING: you are running this with the queryBase parent '
              'class, but only subclasses should be run. Creating a dummy '
              'SQL file.')
        with open(self.config.sql_base + '.sql', 'w') as f:
            f.write('-- Dummy SQL file')

    def submit_job(self, sql_file, output_file='data.dat'):
        '''
        Submits job to the HSC data query service.

        Note that this requires a username and password; these are to be
        provided in the config file. The username should be given directly,
        while the password should first be stored as a environment variable
        and the name of that variable specified in the config file.
        '''
        cf = self.config
        # Basis of the job submission command
        command = f'python hscReleaseQueryDR3.py --user={cf.username} '\
            f'--release-version={cf.release} '\
            f'--password-env={cf.password_env} '\
            f'--format=fits'
        dl_command = ''
        # Use quick mode (short timeout)?
        if cf.do_preview:
            command += ' -p'
        # Download data after querying?
        if cf.do_download:
            if os.path.isfile(output_file):
                print('Found ' + output_file)
                return
            command += ' -d'
            dl_command = ' > ' + output_file
        # Delete submitted job after downloading?
        if cf.do_delete:
            command += ' -D'
        # Add SQL query and output file (if downloading) to command
        command += ' ' + sql_file
        command += dl_command
        # Print and run the command
        print(command)
        os.system(command)

    def run(self):
        '''
        Runs the successive stages of writing and submitting a query.
        '''
        self.write_sql()
        for q, f in zip(self.queries, self.outfiles):
            self.submit_job(q, f)


class queryMetadata(queryBase):
    '''
    Creates queries related to frame metadata.
    '''
    def write_sql(self):
        '''
        Generates queries and writes them to SQL files.
        '''
        cf = self.config
        # Directories for queries and downloaded data
        self.path_queries += 'metadata/'
        self.path_out += 'metadata/'
        # Basis for SQL file names
        sql_base = cf.sql_base
        # Check directories exist
        for p in [self.path_queries, self.path_out]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        self.stout['cols'] = ['*']
        self.stout['from'] = [f'{cf.dr}.frame as frame']
        # Create a query for each field
        for fd in cf.fields:
            ra_min, ra_max, dec_min, dec_max = cf.get_field_boundaries(fd)
            # Account for case where ra_min < ra_max (e.g. SPRING field)
            if ra_max < ra_min:
                stout_ra = f'frame.ra2000 BETWEEN {ra_min} AND 360 OR\n'\
                    f'\tframe.ra2000 BETWEEN 0 AND {ra_max}'
            else:
                stout_ra = f'frame.ra2000 BETWEEN {ra_min} AND {ra_max}'
            stout_dec = f'frame.dec2000 BETWEEN {dec_min} AND {dec_max}'
            self.stout['conds'] = [stout_ra, stout_dec]
            # SQL filename and name for downloaded data
            sql_file_fd = self.path_queries + sql_base + f'_{fd}.sql'
            out_file_fd = self.path_out + sql_base + f'_{fd}.fits'
            # Create queries per band if requested
            if cf.split_by_band:
                self.stout['conds'].append('')  # Dummy to be replaced
                for b in cf.bands.all:
                    stout_b = f'frame.filter=\'{b}\''
                    if b in cf.bands.altnames:
                        for b_alt in cf.bands.altnames[b]:
                            stout_b += f' OR frame.filter=\'{b_alt}\''
                    self.stout['conds'][-1] = stout_b
                    self._assemble_query()
                    # Write to file
                    sql_file_fd_b = sql_file_fd[:-4] + f'_{b}.sql'
                    with open(sql_file_fd_b, 'w') as file:
                        file.write(self.stout_str)
                    # Add to list of queries to submit
                    self.queries.append(sql_file_fd_b)
                    # Add output file name to list
                    out_file_fd_b = f'{out_file_fd[:-4]}_{b}.fits'
                    self.outfiles.append(out_file_fd_b)
            else:
                self._assemble_query()
                # Write to file
                with open(sql_file_fd, 'w') as file:
                    file.write(self.stout_str)
                # Add to list of queries to submit
                self.queries.append(sql_file_fd)
                # Add output file name to list
                self.outfiles.append(out_file_fd)


class queryFlags(queryBase):
    '''
    Creates queries related to quality flags.

    Flags are queried for all sources in the HSC catalogue by default, with
    the option to restrict these to primary detections only via the
    `primary_only` config option. Queries are submitted per field as defined
    in the HSC catalogues.
    '''
    def write_sql(self):
        '''
        Generates queries and writes them to a SQL files.
        '''
        cf = self.config
        # Directories for queries and downloaded data
        self.path_queries += 'flags/'
        self.path_out += 'flags/'
        # Basis for SQL file names
        sql_base = cf.sql_base
        # Check directories exist
        for p in [self.path_queries, self.path_out]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        self.stout['cols'] = [
            'object_id',
            'forced.ra',
            'forced.dec'
        ]
        if not cf.primary_only:
            self.stout['cols'].append('forced.isprimary')

        self.stout['from'] = [f'{cf.dr}.forced as forced']

        for table in cf.flags:
            # Add any other tables to FROM statement
            if table != 'forced':
                self.stout['from'].append(
                    f'{cf.dr}.{table} {table} USING (object_id)'
                )
            # Cycle through requested flags in table
            for col in cf.flags[table]:
                # Cycle through photometric bands for each flag
                for band in cf.flags[table][col]:
                    self.stout['cols'].append(
                        f'{table}.{band}_{col}'
                    )

        # Create a query for each field
        for fd in cf.fields:
            # Get list of subfields belonging to each field
            subs = cf.get_subfields(fd)
            for sfd in subs:
                # Query within current subfield
                self.stout['conds'] = [f'forced.field=\'{sfd}\'']
                if cf.primary_only:
                    self.stout['conds'].append('forced.isprimary=True')

                # Combine all components of query
                self._assemble_query()

                # SQL query file name
                sql_file = f'{self.path_queries}{sql_base}_{fd}_{sfd}.sql'
                # Write to file and append file name to list
                with open(sql_file, 'w') as file:
                    file.write(self.stout_str)
                self.queries.append(sql_file)
                # Output data file name
                out_file = f'{self.path_out}{sql_base}_{fd}_{sfd}.fits'
                self.outfiles.append(out_file)


class queryDustAttenuation(queryBase):
    '''
    Creates queries to retrieve dust attenuation data.

    Data are queried for all sources in the HSC catalogue by default, with
    the option to restrict these to primary detections only via the
    `primary_only` config option. Queries are submitted per field as defined
    in the HSC catalogues.
    '''
    def write_sql(self):
        '''
        Generates queries and writes them to a SQL files.
        '''
        cf = self.config
        # Photometric bands
        bands = cf.bands.all
        # Directories for queries and downloaded data
        self.path_queries += 'dust_attenuation/'
        self.path_out += 'dust_attenuation/'
        # Basis for SQL file names
        sql_base = cf.sql_base
        # Check directories exist
        for p in [self.path_queries, self.path_out]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        self.stout['cols'] = [
            'object_id',
            'forced.ra',
            'forced.dec'
        ] + [
            f'forced.a_{b}' for b in bands
        ]
        if not cf.primary_only:
            self.stout['cols'].append('forced.isprimary')

        self.stout['from'] = [f'{cf.dr}.forced as forced']

        # Create a query for each field
        for fd in cf.fields:
            # Get list of subfields belonging to each field
            subs = cf.get_subfields(fd)
            for sfd in subs:
                # Query within current subfield
                self.stout['conds'] = [
                    f'forced.field=\'{sfd}\'',
                ]
                if cf.primary_only:
                    self.stout['conds'].append('forced.isprimary=True')

                # Combine all components of query
                self._assemble_query()

                # SQL query file name
                sql_file = f'{self.path_queries}{sql_base}_{fd}_{sfd}.sql'
                # Write to file and append file name to list
                with open(sql_file, 'w') as file:
                    file.write(self.stout_str)
                self.queries.append(sql_file)
                # Output data file name
                out_file = f'{self.path_out}{sql_base}_{fd}_{sfd}.fits'
                self.outfiles.append(out_file)


class queryRandoms(queryBase):
    '''
    Creates queries for downloading positions and flags for randoms.

    Positions and the `adjust_density` parameter are queried by default for
    the randoms, with the option to restrict these to primary detections only
    via the `primary_only' config option. Any other columns must be specified
    in the config file - under `band_cols` if they apply to a specific filter,
    or under `extra_cols` otherwise.

    Randoms have a density of 100 per sq. arcmin. Lower densities can be
    requested via the `density` config parameter.
    '''
    def write_sql(self):
        '''
        Generates queries and writes them to a SQL files.
        '''
        cf = self.config
        # Directories for queries and downloaded data
        self.path_queries += 'randoms/'
        self.path_out += 'randoms/'
        # Basis for SQL file names
        sql_base = cf.sql_base
        # Check directories exist
        for p in [self.path_queries, self.path_out]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        self.stout['cols'] = [
            'object_id',
            'ra',
            'dec',
            'adjust_density'
        ]
        if not cf.primary_only:
            self.stout['cols'].append('isprimary')

        # Add each of the requested columns to the list
        if cf.band_cols is not None:
            for col in cf.band_cols:
                for band in cf.band_cols[col]:
                    self.stout['cols'].append(f'{band}_{col}')
        if cf.extra_cols is not None:
            for col in cf.extra_cols:
                self.stout['cols'].append(col)

        self.stout['from'] = [f'{cf.dr}.random']

        # Randoms are downloaded in bins of `adjust_density` with widths of
        # 0.05; figure out how many queries are required to reach the requested
        # density
        ad_w = 0.05
        ad_max = cf.density / 100.
        N_bins = int(ad_max // ad_w) + 1
        ad_bins = [i * ad_w for i in range(N_bins)] + [ad_max]
        # Create a query for each bin
        for i in range(N_bins):
            ad_lo = ad_bins[i]
            ad_hi = ad_bins[i + 1]
            self.stout['conds'] = [
                f'adjust_density >= {ad_lo:.2f}',
                f'adjust_density < {ad_hi:.2f}'
            ]

            # Combine all components of query
            self._assemble_query()

            # SQL query file name
            sql_file = f'{self.path_queries}{sql_base}_{ad_lo:.2f}'\
                f'-{ad_hi:.2f}.sql'
            # Write to file and append file name to list
            with open(sql_file, 'w') as file:
                file.write(self.stout_str)
            self.queries.append(sql_file)
            # Output data file name
            out_file = f'{self.path_out}{sql_base}_{ad_lo:.2f}_-{ad_hi:.2f}'\
                       '.fits'
            self.outfiles.append(out_file)


class queryStarsForDepth(queryBase):
    '''
    Creates queries for downloading data for all stars in the chosen fields.

    This stage will generate and submit queries for downloading data for all
    stars within the chosen fields, with the option to apply additional cuts.
    Additional criteria include:
    - primary detections only;
    - any flags specified in the config must be False;
    - blendedness in the primary band is below a certain threshold.
    By default, will download the RA and Dec., and the flux uncertainties
    corrsponding to the chosen flux type, as these constitute all the info
    required to construct depth maps.
    '''
    def write_sql(self):
        '''
        Generates queries and writes them to a SQL files.
        '''
        cf = self.config
        # Photometric bands
        bands = cf.bands.all
        # Directories for queries and downloaded data
        self.path_queries += 'stars/'
        self.path_out += 'stars/for_depth_map/'
        # Basis for SQL file names
        sql_base = cf.sql_base
        # Check directories exist
        for p in [self.path_queries, self.path_out]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        self.stout['cols'] = [
            'object_id',
            'forced.ra',
            'forced.dec',
        ]
        if not cf.primary_only:
            self.stout['cols'].append('isprimary')

        # For keeping track of which tables need to be joined
        tables = []

        # Fluxes and magnitudes of the chosen type
        tflux = self.flux_tables[cf.mag_type]
        self.stout['cols'].extend(
            [
                f'{tflux}.{b}_{cf.mag_type}_fluxerr'
                for b in bands
            ]
        )
        if tflux != 'forced':
            tables.append(tflux)

        # Conditions for selection (applied to all fields)
        b1 = cf.bands.primary
        stout_cond = [
            f'forced.{b1}_extendedness_value = 0'
        ]
        # Primary sources only?
        if cf.primary_only:
            stout_cond.append(
                'forced.isprimary=True'
            )
        # Flags
        for table in cf.flags:
            for col in cf.flags[table]:
                self.stout['cols'].extend(
                    [f'{table}.{b}_{col}=False'
                        for b in cf.flags[table][col]]
                )
            if table not in tables and table != 'forced':
                tables.append(table)
        # Blendedness cut
        if cf.log_blendedness_max is not None:
            stout_cond.append(
                f'meas2.{b1}_blendedness_abs < '
                f'POWER(10, {cf.log_blendedness_max})'
            )
            tables.append('meas2')

        # Statement specifying the tables to join
        self.stout['from'] = [
            f'{cf.dr}.forced as forced'
        ] + [
            f'{cf.dr}.{table} {table} USING (object_id)'
            for table in tables
        ]

        # Create a query per (sub)field per bin
        for fd in cf.fields:
            # Get list of subfields belonging to each field
            subs = cf.get_subfields(fd)
            for sfd in subs:
                self.stout['conds'] = stout_cond + [f'forced.field=\'{sfd}\'']

                # Combine all components of query
                self._assemble_query()

                # SQL query file name
                sql_file = f'{self.path_queries}{sql_base}_{fd}_{sfd}'\
                           '_fluxerrs.sql'
                # Write to file and append file name to list
                with open(sql_file, 'w') as file:
                    file.write(self.stout_str)
                self.queries.append(sql_file)
                # Output data file name
                out_file = f'{self.path_out}{sql_base}_{fd}_{sfd}'\
                           '_fluxerrs.fits'
                self.outfiles.append(out_file)


class queryGalaxiesBase(queryBase):
    '''
    Base class for querying data for galaxy samples.

    Mostly designed to contain convenience methods specific to querying galaxy
    samples.
    '''
    def __init__(
        self,
        config_file
    ):
        super().__init__(config_file)
        # For keeping track of no. of SQL lines used for sample selection
        self.nlines_samp = 0
        # Paths for additional outputs
        self.path_queries_stars = self.path_queries + 'stars/'
        self.path_out_stars = self.path_out + f'stars/{self.config.run_name}'

    def _write_sql_stars(self):
        '''
        Generates queries for stellar analogues to the target galaxies.

        This creates one query per HSC subfield, and applies all the same cuts
        as for the target galaxy samples with the following exceptions:
        - cuts specific to the galaxy sample are removed (e.g. z cuts for
          for tomographic samples)
        - any extendedness cut is replaced with a cut for point sources
        '''
        cf = self.config
        # Check directories exist
        for p in [self.path_queries_stars, self.path_out_stars]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')
        # Create a query per (sub)field
        for i, fd in enumerate(cf.fields):
            subs = cf.get_subfields(fd)
            for j, sfd in enumerate(subs):
                # Modify the existing conditions list
                if i == j == 0:
                    # Remove any sample-related conditions from the end
                    n = self.nlines_samp
                    if n > 0:
                        self.stout['conds'] = self.stout['conds'][:-n]
                        # Reset nlines_samp to 0
                        self.nlines_samp = 0

                    # Find location of extendedness cut(s)
                    for k, s in enumerate(self.stout['conds']):
                        if 'extendedness_value' in s:
                            self.stout['conds'][k] = \
                                self.stout['conds'][k].replace(
                                    'extendedness_value > 0',
                                    'extendedness_value = 0'
                            )
                # Update the 'field' condition
                self.stout['conds'][-1] = f'forced.field=\'{sfd}\''

                # Combine all components of query
                self._assemble_query()

                # SQL query file name
                sql_file = f'{self.path_queries_stars}stars_{fd}_{sfd}.sql'
                # Write to file and append file name to list
                with open(sql_file, 'w') as file:
                    file.write(self.stout_str)
                self.queries.append(sql_file)
                # Output data file name
                out_file = f'{self.path_out_stars}stars_{fd}_{sfd}.fits'
                self.outfiles.append(out_file)

    def _write_sql_cosmos(self):
        '''
        Generates queries for galaxies in COSMOS field.

        This creates a single query to select galaxies in the COSMOS field
        satisfying all quality control cuts applied to the target samples.
        '''
        # Remove any field and sample conditions from most recent query
        n = self.nlines_samp + 1
        self.stout['conds'] = self.stout['conds'][:-n]
        # Reset nlines_samp to 0
        self.nlines_samp = 0
        # Set field condition to 'cosmos'
        self.stout['conds'].append('forced.field=cosmos')

        # Combine all components of query
        self._assemble_query()

        # SQL query file name
        sql_base = self.config.sql_base
        sql_file = f'{self.path_queries}{sql_base}_cosmos.sql'
        # Write to file and append file name to list
        with open(sql_file, 'w') as file:
            file.write(self.stout_str)
        self.queries.append(sql_file)
        # Output data file name
        out_file = f'{self.path_out}{sql_base}_cosmos.fits'
        self.outfiles.append(out_file)


class queryMaglimTomographic(queryGalaxiesBase):
    '''
    Creates queries for downloading tomographically split, mag-limited samples.

    Given a set of [zmin, zmax] pairs in the config file, this stage will
    generate and submit queries for downloading samples of galaxies with
    photometric redshifts within those bounds. Additional criteria include:
    - primary detections only;
    - any flags specified in the config must be False;
    - blendedness in the primary band is below a certain threshold;
    - dust-corrected primary band magnitude is below a certain threshold;
    - primary band SNR is above a certain threshold;
    - SNR in at least N_sec other bands is above a certain threshold;
    - star/galaxy separation criterion.
    By default, will download the RA and Dec., the fluxes and magnitudes of the
    chosen type (along with their uncertainties), dust attenuation values, the
    best-fit photo-zs and 1-sigma confidence interval for the chosen photo-z
    type, and the corrections required to transform the r- and i-band mags to
    r2 and i2 mags.
    If specified in the config, will also download stars satisfying all of the
    same selection criteria minus the star/galaxy separator.
    '''
    def write_sql(self):
        '''
        Generates queries and writes them to a SQL files.
        '''
        cf = self.config
        # Photometric bands
        bands = cf.bands.all
        # Directories for queries and downloaded data
        self.path_queries += 'galaxies/'
        self.path_out += f'galaxies/{cf.run_name}/'
        paths = [self.path_queries, self.path_out]
        # Basis for SQL file names
        sql_base = cf.sql_base
        # Check directories exist
        for p in paths:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        self.stout['cols'] = [
            'object_id',
            'forced.ra',
            'forced.dec'
        ]
        if not cf.primary_only:
            self.stout['cols'].append('forced.isprimary')
        # For keeping track of which tables need to be joined
        tables = []

        # Fluxes and magnitudes of the chosen type
        tflux = self.flux_tables[cf.mag_type]
        self.stout['cols'].extend(
            [
                f'{tflux}.{b}_{cf.mag_type}_{var}'
                for b in bands
                for var in ['flux', 'fluxerr', 'mag', 'magerr']
            ]
        )
        if tflux != 'forced':
            tables.append(tflux)

        # Dust attenuation values
        self.stout['cols'].extend(
            [f'forced.a_{b}' for b in bands]
        )

        # Photo-z info
        self.stout['cols'].extend(
            [f'{cf.z_table}.photoz_{var}'
             for var in ['best', 'err68_min', 'err68_max']]
        )
        tables.append(cf.z_table)

        # Mag corrections for r and i bands
        self.stout['cols'].extend(
            [f'mag_corr.corr_{b}mag' for b in 'ri']
        )
        tables.append('mag_corr')

        # Extra columns
        if cf.extra_cols is not None:
            for table in cf.extra_cols:
                self.stout['cols'].extend(
                    [f'{table}.{col}' for col in cf.extra_cols[table]]
                )
                if table not in tables and table != 'forced':
                    tables.append(table)

        # Mapping from each band to the corresponding magnitude and flux
        mag_map = {
            b: f'{tflux}.{b}_{cf.mag_type}_mag' for b in bands
        }
        flux_map = {
            b: f'{tflux}.{b}_{cf.mag_type}_flux' for b in bands
        }
        if cf.correct_ri:
            mag_map['r'] += ' - mag_corr.corr_rmag'
            mag_map['i'] += ' - mag_corr.corr_imag'

        # Conditions for selection (applied to all fields and z bins)
        b1 = cf.bands.primary
        b2 = cf.bands.secondary
        stout_cond = [
            f'{mag_map[b1]} - a_{b1} < {cf.maglim}'
        ]
        if cf.primary_only:
            stout_cond.append(
                'forced.isprimary=True'
            )
        # Blendedness cut
        stout_cond.append(
            f'meas2.{b1}_blendedness_abs < '
            f'POWER(10, {cf.log_blendedness_max})'
        )
        tables.append('meas2')
        # SNR cut (primary band)
        stout_cond.append(
            f'{flux_map[b1]} / {flux_map[b1]}err >= {cf.snr_min_primary}',
        )
        # SNR cut (secondary bands)
        stout_cond.append(
            '(' + ' +\n\t '.join(
                [
                    f'(CASE WHEN {flux_map[b]} / {flux_map[b]}err >= '
                    f'{cf.snr_min_secondary} THEN 1 ELSE 0 END)'
                    for b in b2
                ]
            ) + f') >= {cf.N_sec}'
        )
        # Flags
        for table in cf.flags:
            for col in cf.flags[table]:
                stout_cond.extend(
                    [f'{table}.{b}_{col}=False'
                     for b in cf.flags[table][col]]
                )
            if table not in tables and table != 'forced':
                tables.append(table)
        # Star-galaxy separator
        stout_cond.append(
            f'forced.{b1}_extendedness_value > 0'
        )

        # Statement specifying the tables to join
        self.stout['from'] = [
            f'{cf.dr}.forced as forced'
        ] + [
            f'{cf.dr}.{table} {table} USING (object_id)'
            for table in tables
        ]

        # Create a query per (sub)field per bin
        for fd in cf.fields:
            # Get list of subfields belonging to each field
            subs = cf.get_subfields(fd)
            for sfd in subs:
                stout_cond_fd = stout_cond + [f'forced.field=\'{sfd}\'']
                for s in cf.samples:
                    zmin, zmax = cf.samples[s]
                    self.stout['conds'] = stout_cond_fd + [
                        f'{cf.z_table}.photoz_best >= {zmin}',
                        f'{cf.z_table}.photoz_best < {zmax}'
                    ]

                    # Combine all components of query
                    self._assemble_query()

                    # SQL query file name
                    sql_file = f'{self.path_queries}{sql_base}_{fd}_{sfd}_{s}'\
                        '.sql'
                    # Write to file and append file name to list
                    with open(sql_file, 'w') as file:
                        file.write(self.stout_str)
                    self.queries.append(sql_file)
                    # Output data file name
                    out_file = f'{self.path_out}{sql_base}_{fd}_{sfd}_{s}.fits'
                    self.outfiles.append(out_file)
        self.nlines_samp = 2

        # Query galaxies from COSMOS with same quality control applied?
        if cf.query_cosmos:
            self._write_sql_cosmos()

        # Query analogous stars from each (sub)field?
        if cf.query_like_stars:
            self._write_sql_stars()
