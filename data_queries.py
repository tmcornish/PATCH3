##############################################################################
# Classes for querying the HSC database via SQL.
##############################################################################

import os
from stage import baseStage

#############################
#          GLOBALS          #
#############################

# Mapping from different flux types to the tables to which they belong
aper_sizes = [10, 15, 20, 30, 40, 57, 84, 118, 168, 235]
seeings = [0, 1, 2, 3]
conv_sizes = [11, 15, 20, 'kron']
flux_tables = {
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


###############################
#          FUNCTIONS          #
###############################

def _flux_str(b, flux_type):
    '''
    Constructs a string for queries involving fluxes of a specific type.

    Parameters
    ----------
    b: str
        The photometric band of interest.

    flux_type: str
        Type of flux being queried.

    Returns
    -------
    s: str
        String containing the correct table, band and flux type for querying.
    '''
    table = flux_tables[flux_type]
    s = f'{table}.{b}_{flux_type}_flux'
    return s


def _mag_str(b, mag_type, correct_ri=False, correct_dust=False):
    '''
    Constructs a string for queries involving magnitudes of a specific type.

    Parameters
    ----------
    b: str
        The photometric band of interest.

    mag_type: str
        Type of magnitude being queried.

    correct_ri: bool (False)
        Whether or not to apply corrections to the r- and i-band magnitudes
        to acocunt for the mid-survey transition from r/i to r2/i2.

    correct_dust: bool (False)
        Whether or not to include a correction for dust attentuation.

    Returns
        -------
        s: str
            String containing the correct table, band and mag type for
            querying.
    '''
    table = flux_tables[mag_type]
    s = f'{table}.{b}_{mag_type}_mag'
    # r/i correction?
    if correct_ri and b in 'ri':
        s += f' - mag_corr.corr_{b}mag'
    # Dust correction?
    if correct_dust:
        s += f' - forced.a_{b}'
    return s


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
        self.sql = {
            'cols': [],  # Columns to be queried
            'from': [],  # Tables from which they will be queried
            'cond': []   # Conditions to apply when downloading
        }
        self.sql_str = ''
        # Output directory for SQL query files (NOT the data itself)
        self.path_queries = self.config.paths.out + 'sql_queries/'
        # Base directory for the data
        self.path_out = self.config.paths.data
        # Table containing flux/mag info for chosen mag type
        self.tflux = flux_tables[self.config.mag_type]
        # For keeping track of which tables need to be joined
        self.tables_join = []
        # Append with tflux unless 'forced' ('forced' table is always used)
        if self.tflux != 'forced':
            self.tables_join.append(self.tflux)

    def _assemble_query(self):
        '''
        Assembles all components of the query into a single string.
        '''
        sql_cols = 'SELECT ' + ',\n\t'.join(self.sql['cols'])
        sql_from = 'FROM ' + '\n\tLEFT JOIN '.join(self.sql['from'])
        sql_cond = 'WHERE ' + ' AND \n\t'.join(self.sql['cond'])
        self.sql_str = '\n'.join(
            [
                sql_cols,
                sql_from,
                sql_cond
            ]
        ) + '\n;'

    def mag_map(self, b):
        '''
        Constructs a string for queries involving mags of a specific type.

        Parameters
        ----------
        b: str
            The photometric band of interest.

        Returns
        -------
        s: str
            String containing the correct table, band and flux type for
            querying.
        '''
        cf = self.config
        return f'({_mag_str(b, cf.mag_type, cf.correct_ri, cf.correct_dust)})'

    def flux_map(self, b):
        '''
        Constructs a string for queries involving fluxes of a specific type.

        Parameters
        ----------
        b: str
            The photometric band of interest.

        flux_type: str
            Type of flux being queried.

        Returns
        -------
        s: str
            String containing the correct table, band and flux type for
            querying.
        '''
        return _flux_str(b, self.config.mag_type)

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
        self.sql['cols'] = ['*']
        self.sql['from'] = [f'{cf.dr}.frame as frame']
        # Create a query for each field
        for fd in cf.fields:
            ra_min, ra_max, dec_min, dec_max = cf.get_field_boundaries(fd)
            # Account for case where ra_min < ra_max (e.g. SPRING field)
            if ra_max < ra_min:
                sql_ra = f'frame.ra2000 BETWEEN {ra_min} AND 360 OR\n'\
                    f'\tframe.ra2000 BETWEEN 0 AND {ra_max}'
            else:
                sql_ra = f'frame.ra2000 BETWEEN {ra_min} AND {ra_max}'
            sql_dec = f'frame.dec2000 BETWEEN {dec_min} AND {dec_max}'
            self.sql['cond'] = [sql_ra, sql_dec]
            # SQL filename and name for downloaded data
            sql_file_fd = self.path_queries + sql_base + f'_{fd}.sql'
            out_file_fd = self.path_out + sql_base + f'_{fd}.fits'
            # Create queries per band if requested
            if cf.split_by_band:
                self.sql['cond'].append('')  # Dummy to be replaced
                for b in cf.bands.all:
                    sql_b = f'frame.filter=\'{b}\''
                    if b in cf.bands.altnames:
                        for b_alt in cf.bands.altnames[b]:
                            sql_b += f' OR frame.filter=\'{b_alt}\''
                    self.sql['cond'][-1] = sql_b
                    self._assemble_query()
                    # Write to file
                    sql_file_fd_b = sql_file_fd[:-4] + f'_{b}.sql'
                    with open(sql_file_fd_b, 'w') as file:
                        file.write(self.sql_str)
                    # Add to list of queries to submit
                    self.queries.append(sql_file_fd_b)
                    # Add output file name to list
                    out_file_fd_b = f'{out_file_fd[:-5]}_{b}.fits'
                    self.outfiles.append(out_file_fd_b)
            else:
                self._assemble_query()
                # Write to file
                with open(sql_file_fd, 'w') as file:
                    file.write(self.sql_str)
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
        self.path_out += f'flags/{cf.run_name}/'
        # Basis for SQL file names
        sql_base = cf.sql_base
        # Check directories exist
        for p in [self.path_queries, self.path_out]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        self.sql['cols'] = [
            'object_id',
            'forced.ra',
            'forced.dec'
        ]
        if not cf.primary_only:
            self.sql['cols'].append('forced.isprimary')

        self.sql['from'] = [f'{cf.dr}.forced as forced']

        for table in cf.flags:
            # Add any other tables to FROM statement
            if table != 'forced':
                self.sql['from'].append(
                    f'{cf.dr}.{table} {table} USING (object_id)'
                )
            # Cycle through requested flags in table
            for col in cf.flags[table]:
                # Cycle through photometric bands for each flag
                for band in cf.flags[table][col]:
                    self.sql['cols'].append(
                        f'{table}.{band}_{col}'
                    )

        # Create a query for each field
        for fd in cf.fields:
            # Get list of subfields belonging to each field
            subs = cf.get_subfields(fd)
            for sfd in subs:
                # Query within current subfield
                self.sql['cond'] = [f'forced.field=\'{sfd}\'']
                if cf.primary_only:
                    self.sql['cond'].append('forced.isprimary=True')

                # Combine all components of query
                self._assemble_query()

                # SQL query file name
                sql_file = f'{self.path_queries}{sql_base}_{fd}_{sfd}.sql'
                # Write to file and append file name to list
                with open(sql_file, 'w') as file:
                    file.write(self.sql_str)
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
        self.sql['cols'] = [
            'object_id',
            'forced.ra',
            'forced.dec'
        ] + [
            f'forced.a_{b}' for b in bands
        ]
        if not cf.primary_only:
            self.sql['cols'].append('forced.isprimary')

        self.sql['from'] = [f'{cf.dr}.forced as forced']

        # Create a query for each field
        for fd in cf.fields:
            # Get list of subfields belonging to each field
            subs = cf.get_subfields(fd)
            for sfd in subs:
                # Query within current subfield
                self.sql['cond'] = [
                    f'forced.field=\'{sfd}\'',
                ]
                if cf.primary_only:
                    self.sql['cond'].append('forced.isprimary=True')

                # Combine all components of query
                self._assemble_query()

                # SQL query file name
                sql_file = f'{self.path_queries}{sql_base}_{fd}_{sfd}.sql'
                # Write to file and append file name to list
                with open(sql_file, 'w') as file:
                    file.write(self.sql_str)
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
        self.sql['cols'] = [
            'object_id',
            'ra',
            'dec',
            'adjust_density'
        ]
        if not cf.primary_only:
            self.sql['cols'].append('isprimary')

        # Add each of the requested columns to the list
        if cf.band_cols is not None:
            for col in cf.band_cols:
                for band in cf.band_cols[col]:
                    self.sql['cols'].append(f'{band}_{col}')
        if cf.extra_cols is not None:
            for col in cf.extra_cols:
                self.sql['cols'].append(col)

        self.sql['from'] = [f'{cf.dr}.random']

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
            self.sql['cond'] = [
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
                file.write(self.sql_str)
            self.queries.append(sql_file)
            # Output data file name
            out_file = f'{self.path_out}{sql_base}_{ad_lo:.2f}-{ad_hi:.2f}'\
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
        self.sql['cols'] = [
            'object_id',
            'forced.ra',
            'forced.dec',
        ]
        if not cf.primary_only:
            self.sql['cols'].append('isprimary')

        # Fluxes and magnitudes of the chosen type
        self.sql['cols'].extend(
            [
                f'{self.flux_map(b)}err'
                for b in bands
            ]
        )

        # Conditions for selection (applied to all fields)
        b1 = cf.bands.primary
        sql_cond = [
            f'forced.{b1}_extendedness_value = 0'
        ]
        # Primary sources only?
        if cf.primary_only:
            sql_cond.append(
                'forced.isprimary=True'
            )
        # Flags
        for table in cf.flags:
            for col in cf.flags[table]:
                self.sql['cols'].extend(
                    [f'{table}.{b}_{col}=False'
                        for b in cf.flags[table][col]]
                )
            if table not in self.tables_join and table != 'forced':
                self.tables_join.append(table)
        # Blendedness cut
        if cf.log_blendedness_max is not None:
            sql_cond.append(
                f'meas2.{b1}_blendedness_abs < '
                f'POWER(10, {cf.log_blendedness_max})'
            )
            self.tables_join.append('meas2')

        # Statement specifying the tables to join
        self.sql['from'] = [
            f'{cf.dr}.forced as forced'
        ] + [
            f'{cf.dr}.{table} {table} USING (object_id)'
            for table in self.tables_join
        ]

        # Create a query per (sub)field per bin
        for fd in cf.fields:
            # Get list of subfields belonging to each field
            subs = cf.get_subfields(fd)
            for sfd in subs:
                self.sql['cond'] = sql_cond + [f'forced.field=\'{sfd}\'']

                # Combine all components of query
                self._assemble_query()

                # SQL query file name
                sql_file = f'{self.path_queries}{sql_base}_{fd}_{sfd}'\
                           '_fluxerrs.sql'
                # Write to file and append file name to list
                with open(sql_file, 'w') as file:
                    file.write(self.sql_str)
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

    def _parse_cuts(self, sample, sql_cond):
        '''
        Constructs SQL-formatted cuts from the config file.

        By default, handles True and False conditions. Additional cuts can
        be specified in the config for subclasses but their handling must
        also be defined for each subclass separately.

        Parameters
        ----------
        sample: str
            Name of the sample for which the cuts are intended.
            Can be 'all' if cuts are to be applied regardless of sample.

        sql_cond: list[str]
            List of conditions to which the new cuts will be appended.
        '''
        cf = self.config
        cond = []
        # Columns whose values should be False
        for t in cf.cuts.false_cols[sample]:
            for col in cf.cuts.false_cols[sample][t]:
                if col == 'merge_peak':
                    cond.extend(
                        [f'{t}.{col}_{b}=False'
                            for b in cf.cuts.false_cols[sample][t][col]]
                    )
                else:
                    cond.extend(
                        [f'{t}.{b}_{col}=False'
                            for b in cf.cuts.false_cols[sample][t][col]]
                    )
            if t not in self.tables_join and t != 'forced':
                self.tables_join.append(t)
        # Columns whose values should be True
        for t in cf.cuts.true_cols[sample]:
            for col in cf.cuts.true_cols[sample][t]:
                if col == 'merge_peak':
                    cond.extend(
                        [f'{t}.{col}_{b}=True'
                            for b in cf.cuts.true_cols[sample][t][col]]
                    )
                else:
                    cond.extend(
                        [f'{t}.{b}_{col}=True'
                            for b in cf.cuts.true_cols[sample][t][col]]
                    )
            if t not in self.tables_join and t != 'forced':
                self.tables_join.append(t)
        sql_cond.extend(cond)

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
                        self.sql['cond'] = self.sql['cond'][:-n]
                        # Reset nlines_samp to 0
                        self.nlines_samp = 0

                    # Find location of extendedness cut(s)
                    for k, s in enumerate(self.sql['cond']):
                        if 'extendedness_value' in s:
                            self.sql['cond'][k] = \
                                self.sql['cond'][k].replace(
                                    'extendedness_value > 0',
                                    'extendedness_value = 0'
                            )
                # Update the 'field' condition
                self.sql['cond'][-1] = f'forced.field=\'{sfd}\''

                # Combine all components of query
                self._assemble_query()

                # SQL query file name
                sql_file = f'{self.path_queries_stars}stars_{fd}_{sfd}.sql'
                # Write to file and append file name to list
                with open(sql_file, 'w') as file:
                    file.write(self.sql_str)
                self.queries.append(sql_file)
                # Output data file name
                out_file = f'{self.path_out_stars}stars_{fd}_{sfd}.fits'
                self.outfiles.append(out_file)

    def _write_sql_cosmos_dud(self, sample):
        '''
        Generates queries for galaxies in COSMOS deep/ultra-deep catalogue.

        This creates a single query to select galaxies in the COSMOS field
        satisfying all quality control and sample selection cuts applied
        to the target samples.
        '''
        # Replace the field condition
        n = self.nlines_samp + 1
        self.sql['cond'][-n] = 'forced.field=\'cosmos\''

        # Combine all components of query
        self._assemble_query()
        # Replace all instances of the Wide data release with 'pdr3_dud'
        self.sql_str = self.sql_str.replace(
            self.config.dr,
            self.config.dr_cosmos
        )

        # SQL query file name
        sql_base = self.config.sql_base
        sql_file = f'{self.path_queries}{sql_base}_cosmos_{sample}.sql'
        # Write to file and append file name to list
        with open(sql_file, 'w') as file:
            file.write(self.sql_str)
        self.queries.append(sql_file)
        # Output data file name
        out_file = f'{self.path_out}{sql_base}_cosmos_{sample}.fits'
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
        # Number of lines specific to each sample
        self.nlines_samp = 2
        # Check directories exist
        for p in paths:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        self.sql['cols'] = [
            'object_id',
            'forced.ra',
            'forced.dec'
        ]
        if not cf.primary_only:
            self.sql['cols'].append('forced.isprimary')

        # Fluxes and magnitudes of the chosen type
        self.sql['cols'].extend(
            [
                f'{self.tflux}.{b}_{cf.mag_type}_{var}'
                for b in bands
                for var in ['flux', 'fluxerr', 'mag', 'magerr']
            ]
        )

        # Dust attenuation values
        self.sql['cols'].extend(
            [f'forced.a_{b}' for b in bands]
        )

        # Photo-z info
        self.sql['cols'].extend(
            [f'{cf.z_table}.photoz_{var}'
             for var in ['best', 'err68_min', 'err68_max']]
        )
        self.tables_join.append(cf.z_table)

        # Mag corrections for r and i bands
        self.sql['cols'].extend(
            [f'mag_corr.corr_{b}mag' for b in 'ri']
        )
        self.tables_join.append('mag_corr')

        # Extra columns
        if cf.extra_cols is not None:
            for table in cf.extra_cols:
                self.sql['cols'].extend(
                    [f'{table}.{col}' for col in cf.extra_cols[table]]
                )
                if table not in self.tables_join and table != 'forced':
                    self.tables_join.append(table)

        # Conditions for selection (applied to all fields and z bins)
        b1 = cf.bands.primary
        b2 = cf.bands.secondary
        sql_cond = [
            self.mag_map(b1) + f' < {cf.maglim}'
        ]
        if cf.primary_only:
            sql_cond.append(
                'forced.isprimary=True'
            )
        # Blendedness cut
        sql_cond.append(
            f'meas2.{b1}_blendedness_abs < '
            f'POWER(10, {cf.log_blendedness_max})'
        )
        self.tables_join.append('meas2')
        # SNR cut (primary band)
        sql_cond.append(
            f'{self.flux_map(b1)} / {self.flux_map(b1)}err >= '
            f'{cf.snr_min_primary}',
        )
        # SNR cut (secondary bands)
        sql_cond.append(
            '(' + ' +\n\t '.join(
                [
                    f'(CASE WHEN {self.flux_map(b)} / {self.flux_map(b)}err '
                    f'>= {cf.snr_min_secondary} THEN 1 ELSE 0 END)'
                    for b in b2
                ]
            ) + f') >= {cf.N_sec}'
        )
        # Flags
        self._parse_cuts('all', sql_cond)
        # Star-galaxy separator
        sql_cond.append(
            f'forced.{b1}_extendedness_value > 0'
        )

        # Statement specifying the tables to join
        self.sql['from'] = [
            f'{cf.dr}.forced as forced'
        ] + [
            f'{cf.dr}.{table} {table} USING (object_id)'
            for table in self.tables_join
        ]

        # Create a query per (sub)field per bin
        for i_fd, fd in enumerate(cf.fields):
            # Get list of subfields belonging to each field
            subs = cf.get_subfields(fd)
            for sfd in subs:
                sql_cond_fd = sql_cond + [f'forced.field=\'{sfd}\'']
                for s in cf.samples:
                    zmin, zmax = cf.samples[s]
                    self.sql['cond'] = sql_cond_fd + [
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
                        file.write(self.sql_str)
                    self.queries.append(sql_file)
                    # Output data file name
                    out_file = f'{self.path_out}{sql_base}_{fd}_{sfd}_{s}.fits'
                    self.outfiles.append(out_file)

                    # Query COSMOS galaxies with same cuts applied?
                    if cf.query_cosmos and i_fd == 0:
                        self._write_sql_cosmos_dud(s)

        # Query analogous stars from each (sub)field?
        if cf.query_like_stars:
            self._write_sql_stars()


class queryLBGs(queryGalaxiesBase):
    '''
    Creates queries for downloading various LBG dropout samples.

    Applies the colour-colour and colour-magnitude cuts employed in the
    GOLDRUSH analyses (Harikane et al. 2022) to select g-, r-, i-, and/or
    z-dropouts. Flags can also be specified in the config file for each
    sample.
    '''
    def __init__(
            self,
            config_file
    ):
        super().__init__(config_file)
        # Detection bands for each dropout sample
        self.det_bands = {
            'g_dropouts': 'i',
            'r_dropouts': 'z',
            'i_dropouts': 'z',
            'z_dropouts': 'y',
        }

    def mag_map_snr(self, b):
        '''
        Similar to `mag_map()`, but using the SNR magnitude type.

        Parameters
        ----------
        b: str
            The photometric band of interest.

        Returns
        -------
        s: str
            String containing the correct table, band and flux type for
            querying.
        '''
        cf = self.config
        s = f'({_mag_str(b, cf.mag_type_snr, cf.correct_ri, cf.correct_dust)})'
        return s

    def flux_map_snr(self, b):
        '''
        Similar to `flux_map()`, but using the SNR magnitude type.

        Parameters
        ----------
        b: str
            The photometric band of interest.

        Returns
        -------
        s: str
            String containing the correct table, band and flux type for
            querying.
        '''
        return _flux_str(b, self.config.mag_type_snr)

    def _parse_cuts(self, sample, sql_cond):
        '''
        Constructs SQL-formatted cuts from the config file.

        Parameters
        ----------
        sample: str
            Name of the sample for which the cuts are intended.
            Must be one of ['all', 'g_dropouts', 'r_dropouts',
            'i_dropouts', 'z_dropouts'].

        sql_cond: list[str]
            List of conditions to which the new cuts will be appended.
        '''
        super()._parse_cuts(sample, sql_cond)

        cf = self.config
        cond = []
        if sample != 'all':
            # Blendedness cuts
            t = 'meas2'
            cond.extend(
                [f'{t}.{b}_blendedness_abs < {cf.blendedness_abs_max}'
                 for b in cf.cuts.blendedness[sample]]
            )
            if t not in self.tables_join:
                self.tables_join.append(t)
        sql_cond.extend(cond)

    def _dropout_conditions(self, sample, sql_cond):
        '''
        Constructs the colour-colour/magnitude cuts for a dropout sample.

        Parameters
        ----------
        sample: str
            Dropout band ('g_dropouts', 'r_dropouts', 'i_dropouts', or
            'z_dropouts').

        sql_cond: list[str]
            List of SQL-formatted conditions.
        '''
        # g-dropouts
        if sample == 'g_dropouts':
            cond = [
                f'({self.flux_map_snr("i")} / {self.flux_map_snr("i")}err) '
                '>= 5',
                f'({self.mag_map("g")} - {self.mag_map("r")}) > 1.0',
                f'({self.mag_map("r")} - {self.mag_map("i")}) < 1.0',
                f'({self.mag_map("g")} - {self.mag_map("r")}) > 1.5 '
                f'* ({self.mag_map("r")} - {self.mag_map("i")}) + 0.8'
            ]
        # r-dropouts
        elif sample == 'r_dropouts':
            cond = [
                f'({self.flux_map_snr("g")} / {self.flux_map_snr("g")}err) '
                '< 2.0',
                f'({self.flux_map_snr("z")} / {self.flux_map_snr("z")}err) '
                '>= 5',
                f'({self.mag_map("r")} - {self.mag_map("i")}) > 1.2',
                f'({self.mag_map("i")} - {self.mag_map("z")}) < 0.7',
                f'({self.mag_map("r")} - {self.mag_map("i")}) > 1.5 '
                f'* ({self.mag_map("i")} - {self.mag_map("z")}) + 1.0'
            ]
        # i-dropouts
        elif sample == 'i_dropouts':
            cond = [
                f'({self.flux_map_snr("g")} / {self.flux_map_snr("g")}err) '
                '< 2.0',
                f'({self.flux_map_snr("r")} / {self.flux_map_snr("r")}err) '
                '< 2.0',
                f'({self.flux_map_snr("z")} / {self.flux_map_snr("z")}err) '
                '>= 5',
                f'({self.flux_map_snr("y")} / {self.flux_map_snr("y")}err) '
                '>= 4',
                f'({self.mag_map("i")} - {self.mag_map("z")}) > 1.5',
                f'({self.mag_map("z")} - {self.mag_map("y")}) < 0.7',
                f'({self.mag_map("i")} - {self.mag_map("z")}) > 2.0 '
                f'* ({self.mag_map("z")} - {self.mag_map("y")}) + 1.1'
            ]
        # z-dropouts
        elif sample == 'z_dropouts':
            cond = [
                f'({self.flux_map_snr("g")} / {self.flux_map_snr("g")}err) '
                '< 2.0',
                f'({self.flux_map_snr("r")} / {self.flux_map_snr("r")}err) '
                '< 2.0',
                f'({self.flux_map_snr("i")} / {self.flux_map_snr("i")}err) '
                '< 2.0',
                f'({self.flux_map_snr("y")} / {self.flux_map_snr("y")}err) '
                '>= 4',
                f'({self.mag_map("z")} - {self.mag_map("y")}) > 1.6'
            ]
        # Limiting magnitude in detection band
        mlim = self.config.cuts.maglim[sample]
        if mlim is not None:
            cond.append(
                f'{self.mag_map_snr(self.det_bands[sample])} < {mlim}'
            )
        sql_cond.extend(cond)

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
        self.sql['cols'] = [
            'object_id',
            'forced.ra',
            'forced.dec'
        ]
        if not cf.primary_only:
            self.sql['cols'].append('forced.isprimary')

        # Magnitudes of the chosen type for colour cuts
        self.sql['cols'].extend(
            [
                f'{_mag_str(b, cf.mag_type)}{var}'
                for b in bands
                for var in ['', 'err']
            ]
        )
        # Fluxes of the chosen type for SNR cuts
        self.sql['cols'].extend(
            [
                f'{_flux_str(b, cf.mag_type_snr)}{var}'
                for b in bands
                for var in ['', 'err']
            ]
        )
        tflux_snr = flux_tables[cf.mag_type_snr]
        if tflux_snr != 'forced':
            self.tables_join.append(tflux_snr)

        # Dust attenuation values
        self.sql['cols'].extend(
            [f'forced.a_{b}' for b in bands]
        )

        # Photo-z info
        self.sql['cols'].extend(
            [f'{cf.z_table}.photoz_{var}'
                for var in ['best', 'err68_min', 'err68_max']]
        )
        self.tables_join.append(cf.z_table)

        # Mag corrections for r and i bands
        self.sql['cols'].extend(
            [f'mag_corr.corr_{b}mag' for b in 'ri']
        )
        self.tables_join.append('mag_corr')

        # Extra columns
        if cf.extra_cols is not None:
            for table in cf.extra_cols:
                self.sql['cols'].extend(
                    [f'{table}.{col}' for col in cf.extra_cols[table]]
                )
                if table not in self.tables_join and table != 'forced':
                    self.tables_join.append(table)

        # True/False and blendedness cuts applied to all samples
        sql_cond = []
        self._parse_cuts('all', sql_cond)
        # Cuts on inputcount_value in each band, applied to all samples
        sql_cond.extend(
            [f'forced.{b}_inputcount_value >= {cf.inputcount_min[b]}'
             for b in bands]
        )

        # Create a query per (sub)field per sample
        for i_fd, fd in enumerate(cf.fields):
            # Get list of subfields belonging to each field
            subs = cf.get_subfields(fd)
            for sfd in subs:
                sql_cond_fd = sql_cond + [f'forced.field=\'{sfd}\'']
                # Cuts specific to each dropout sample
                for b in 'griz':
                    sample = f'{b}_dropouts'
                    sql_cond_samp = sql_cond_fd.copy()
                    self._parse_cuts(sample, sql_cond_samp)
                    self._dropout_conditions(sample, sql_cond_samp)

                    self.sql['cond'] = sql_cond_samp
                    self.nlines_samp = len(self.sql['cond']) - len(sql_cond_fd)

                    # Statement specifying the tables to join
                    self.sql['from'] = [
                        f'{cf.dr}.forced as forced'
                    ] + [
                        f'{cf.dr}.{table} {table} USING (object_id)'
                        for table in self.tables_join
                    ]

                    # Combine all components of query
                    self._assemble_query()

                    # SQL query file name
                    sql_file = f'{self.path_queries}{sql_base}_{fd}_{sfd}'\
                        f'_{sample}.sql'
                    # Write to file and append file name to list
                    with open(sql_file, 'w') as file:
                        file.write(self.sql_str)
                    self.queries.append(sql_file)
                    # Output data file name
                    out_file = f'{self.path_out}{sql_base}_{fd}_{sfd}'\
                        f'_{sample}.fits'
                    self.outfiles.append(out_file)

                    # Query COSMOS galaxies with same cuts?
                    if cf.query_cosmos and i_fd == 0:
                        self._write_sql_cosmos_dud(sample)

        # Query analogous stars from each (sub)field?
        if cf.query_like_stars:
            self._write_sql_stars()
