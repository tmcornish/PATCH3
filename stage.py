##############################################################################
# Basic class from which all stages of the pipeline inherit.
##############################################################################

import os
from configuration import PipelineConfig as PC


class baseStage(object):
    '''
    Base stage upon which all others are built.

    Parameters
    ----------
    config_file: str
        Path to a pipeline YAML config file.
    '''
    def __init__(
        self,
        config_file
    ):
        self.config = PC(config_file, stage=self.__class__.__name__)
        self.config_file = config_file

    def run(self):
        '''
        Contains all code that will be run upon execution of the stage.

        Here this is just a placeholder, as it will be defined individually
        for any children classes.
        '''
        return

    def execute(self):
        '''
        Method for executing the stage.

        This method depends on the machine on which the pipeline is being run.
        For example, if being run locally, the `run()` method is simply run
        as-is, but if being run on the Imperial HPC, this will generate qsub
        job scripts and submit those. NOTE: So far this only has support for
        execution locally or on qsub-based clusters like the Imperial HPC.
        '''
        cf = self.config
        platform = cf.platform
        stagename = self.__class__.__name__

        if platform == 'local':
            import cusp.output_utils as ou
            print(
                ou.colour_string(
                    ou.string_important(stagename),
                    'orange'
                ) + '\n'
            )
            self.run()

        elif platform == 'imperial-cx3':
            print(f'Submitting job for stage {stagename}...')
            # Retrieve relevant information from config
            walltime = cf.qsub.walltime
            resources = cf.qsub.resources
            conda_env = cf.conda_env
            path_pipe = cf.paths.pipeline
            runfile = cf.runfile
            envsfile = cf.envsfile
            # Assemble a bash script for submitting this stage as a job
            res_str = ':'.join([f'{k}={resources[k]}' for k in resources])
            py_str = [
                f'from {runfile[:-3]} import {stagename}',
                f's = {stagename}(\'{self.config_file}\')',
                's.run()'
            ]
            py_str = '; '.join(py_str)
            script = [
                '#!/bin/bash',
                f'#PBS -l {res_str}',
                f'#PBS -l walltime={walltime}',
                f'#PBS -N PATCH3-{stagename}',
                '',
                'eval "$(~/miniforge3/bin/conda shell.bash hook)"',
                f'conda activate {conda_env}',
                '',
                f'cd {path_pipe}',
                '',
                f'source {envsfile}',
                f'python -c "{py_str}"'
                ]
            script = '\n'.join(script)

            # Write script to a file
            path_scripts = path_pipe + 'qsub_scripts/'
            script_file = path_scripts + f'{stagename}.sh'
            if not os.path.exists(path_scripts):
                os.mkdir(path_scripts)
            with open(script_file, 'w') as f:
                f.write(script)

            # Check if previous job has been submitted
            depend = ''
            prevjob_file = path_pipe + f'prevjob_{cf.run_name}.txt'
            if os.path.exists(prevjob_file):
                with open(prevjob_file, 'r') as f:
                    jobid = f.read().strip('\n')
                    depend = f'-W depend=afterok:{jobid} '
            # Submit job to queue
            os.system(f'qsub {depend}{script_file} > {prevjob_file}')

        else:
            raise ValueError(
                'PATCH3 is currently only supported for running locally or on '
                'the Imperial College London HPC. Support for SLURM-based '
                'clusters like NERSC will be added in future updates.'
            )
