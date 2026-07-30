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
        baseStage.__init__(self, config_file)
        self.queries = []
        self.outfiles = []
        # Output directory fro SQL query files (NOT the data itself)
        self.path_queries = self.config.paths.out + 'sql_queries/'
        # Length of extension for chosen output format
        self.n_ext = len(self.config.format) + 1

    def write_sql(self):
        '''
        Generates a query and writes it to an SQL file.

        Here this is just a placeholder; it is properly defined for
        each subclass.
        '''
        print('WARNING: you are running this with the queryBase parent '
              'class, but only subclasses should be run. Creating a dummy '
              'SQL file.')
        with open(self.config.sql_base, 'w') as f:
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
            f'--format={cf.format}'
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
        path_queries = self.path_queries + 'metadata/'
        path_out = cf.paths.data + 'metadata/'
        # Basis for SQL file names
        sql_base = cf.sql_base
        # Check directories exist
        for p in [path_queries, path_out]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        stout_base = [
            'SELECT *',
            f'FROM {cf.dr}.frame as frame',
            'WHERE '
        ]
        stout_base = '\n'.join(stout_base)
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
            stout_fd = [
                stout_ra,
                stout_dec
            ]
            stout_fd = ' AND \n\t'.join(stout_fd)
            stout_fd = stout_base + stout_fd
            # SQL filename and name for downloaded data
            sql_file_fd = path_queries + sql_base[:-4] + f'_{fd}.sql'
            out_file_fd = path_out + sql_base[:-4] + f'_{fd}.{cf.format}'
            # Create queries per band if requested
            if cf.split_by_band:
                for b in cf.bands.all:
                    stout_b = f' AND \n\tframe.filter=\'{b}\''
                    if b in cf.bands.altnames:
                        for b_alt in cf.bands.altnames[b]:
                            stout_b += f' OR frame.filter=\'{b_alt}\''
                    stout = stout_fd + stout_b + '\n;'
                    # Write to file
                    sql_file_fd_b = sql_file_fd[:-4] + f'_{b}.sql'
                    with open(sql_file_fd_b, 'w') as file:
                        file.write(stout)
                    # Add to list of queries to submit
                    self.queries.append(sql_file_fd_b)
                    # Add output file name to list
                    out_file_fd_b = f'{out_file_fd[:-self.n_ext]}'\
                        f'_{b}.{cf.format}'
                    self.outfiles.append(out_file_fd_b)
            else:
                stout = stout_fd + '\n;'
                # Write to file
                with open(sql_file_fd, 'w') as file:
                    file.write(stout)
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
        path_queries = self.path_queries + 'flags/'
        path_out = cf.paths.data + 'flags/'
        # Basis for SQL file names
        sql_base = cf.sql_base
        # Check directories exist
        for p in [path_queries, path_out]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        stout_cols = [
            'SELECT object_id',
            'forced.ra',
            'forced.dec'
        ]
        if not cf.primary_only:
            stout_cols.append('forced.isprimary')

        stout_from = [
            f'FROM {cf.dr}.forced as forced'
        ]

        for table in cf.flags:
            # Add any other tables to FROM statement
            if table != 'forced':
                stout_from.append(
                    f'{cf.dr}.{table} {table} USING (object_id)'
                )
            # Cycle through requested flags in table
            for col in cf.flags[table]:
                # Cycle through photometric bands for each flag
                for band in cf.flags[table][col]:
                    stout_cols.append(
                        f'{table}.{band}_{col}'
                    )
        stout_cols = ',\n\t'.join(stout_cols)
        stout_from = '\n\tLEFT JOIN '.join(stout_from)

        # Create a query for each field
        for fd in cf.fields:
            # Get list of subfields belonging to each field
            subs = cf.get_subfields(fd)
            for sfd in subs:
                # Query within current subfield
                stout_cond = [
                    f'WHERE forced.field=\'{sfd}\'',
                ]
                if cf.primary_only:
                    stout_cond.append('forced.isprimary=True')
                stout_cond = ' AND \n\t'.join(stout_cond)

                # Combine all components of query
                stout = [
                    stout_cols,
                    stout_from,
                    stout_cond
                ]
                stout = '\n'.join(stout) + '\n;'

                # SQL query file name
                sql_file = f'{path_queries}{sql_base[:-4]}_{fd}_{sfd}.sql'
                # Write to file and append file name to list
                with open(sql_file, 'w') as file:
                    file.write(stout)
                self.queries.append(sql_file)
                # Output data file name
                out_file = f'{path_out}{sql_base[:-4]}_{fd}_{sfd}.{cf.format}'
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
        path_queries = self.path_queries + 'randoms/'
        path_out = cf.paths.data + 'randoms/'
        # Basis for SQL file names
        sql_base = cf.sql_base
        # Check directories exist
        for p in [path_queries, path_out]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Begin assembling query
        stout_cols = [
            'SELECT object_id',
            'ra',
            'dec',
            'adjust_density'
        ]
        if not cf.primary_only:
            stout_cols.append('isprimary')

        # Add each of the requested columns to the list
        if cf.band_cols is not None:
            for col in cf.band_cols:
                for band in cf.band_cols[col]:
                    stout_cols.append(f'{band}_{col}')
        if cf.extra_cols is not None:
            for col in cf.extra_cols:
                stout_cols.append(col)
        stout_cols = ',\n\t'.join(stout_cols)

        stout_from = f'FROM {cf.dr}.random'

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
            stout_cond = (
                f'WHERE adjust_density >= {ad_lo} AND \n'
                f'\tadjust_density < {ad_hi}'
            )

            # Combine all components of query
            stout = [
                stout_cols,
                stout_from,
                stout_cond
            ]
            stout = '\n'.join(stout) + '\n;'

            # SQL query file name
            sql_file = f'{path_queries}{sql_base[:-4]}_{ad_lo:.2f}'\
                f'-{ad_hi:.2f}.sql'
            # Write to file and append file name to list
            with open(sql_file, 'w') as file:
                file.write(stout)
            self.queries.append(sql_file)
            # Output data file name
            out_file = f'{path_out}{sql_base[:-4]}_{ad_lo:.2f}_-{ad_hi:.2f}'\
                f'.{cf.format}'
            self.outfiles.append(out_file)
