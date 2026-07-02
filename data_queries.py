##############################################################################
# Classes for querying the HSC database via SQL.
##############################################################################

import os
from configuration import PipelineConfig as PC


############################
#          STAGES          #
############################

class queryBase(object):
    '''
    Base class for generating queries; should not be used directly.
    '''
    def __init__(
        self,
        config_file,
    ):
        cf = PC(config_file, stage=self.__class__.__name__)
        self.config = cf
        self.queries = []
        self.outfiles = []
        self.path_queries = cf.paths.out + 'sql_queries/'

    def write_sql(self, sql_file):
        '''
        Generates a query and writes it to an SQL file.

        Here this is just a placeholder; it is properly defined for
        each subclass.
        '''
        print('WARNING: you are running this with the queryBase parent '
              'class, but only subclasses should be run. Creating a dummy '
              'SQL file.')
        with open(sql_file, 'w') as f:
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
            f'--release-version={cf.release} --password-env={cf.password_env}'\
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
    def write_sql(self, sql_file):
        '''
        Generates a query and writes it to an SQL file.
        '''
        cf = self.config
        # Directories for queries and downloaded data
        path_queries = self.path_queries + 'metadata/'
        path_out = cf.paths.data + 'metadata/'
        # Check directories exist
        for p in [path_queries, path_out]:
            if not os.path.exists(p):
                os.system(f'mkdir -p {p}')

        # Length of extension for chosen output format
        n_ext = len(cf.format) + 1

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
            sql_file_fd = path_queries + sql_file[:-4] + f'_{fd}.sql'
            out_file_fd = path_out + sql_file[:-4] + f'_{fd}.{cf.format}'
            # Create queries per band if requested
            if cf.split_by_band:
                for b in cf.bands.all:
                    stout_b = f' AND \n\tframe.band=\'{b}\''
                    stout = stout_fd + stout_b + '\n;'
                    # Write to file
                    sql_file_fd_b = sql_file_fd[:-4] + f'_{b}.sql'
                    with open(sql_file_fd_b, 'w') as file:
                        file.write(stout)
                    # Add to list of queries to submit
                    self.queries.append(sql_file_fd_b)
                    # Add output file name to list
                    out_file_fd_b = out_file_fd[:-n_ext] + f'_{b}.{cf.format}'
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
