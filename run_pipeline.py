###########################################################################
# Master script for running different stages of the pipeline.
###########################################################################

# Import necessary packages/modules
import sys
import os
import yaml

# Retrieve the path of the pipeline config file
config_file = sys.argv[1]

# Identify which stages need to be run
with open(config_file) as f:
    config = yaml.safe_load(f)
stages = config['stages']

# Cycle through stages and run if told to do so
for stage in stages:
    if stages[stage]:
        config_stage = config[stage]
        runfile = config_stage['runfile']

        # Run execute() method of the current stage
        py_str = [
            f'from {runfile[:-3]} import {stage}',
            f's = {stage}(\'{config_file}\')',
            's.execute()'
        ]
        py_str = '; '.join(py_str)
        os.system(f'python -c "{py_str}"')
