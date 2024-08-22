#!/usr/bin/env python3

import pandas as pd
import os
import yaml

# Get output folder using snakemake or, if running script independently, directly from config file
try:
    output = snakemake.params.output
except NameError:
    with open('config/config.yaml', 'r') as f:
        config = yaml.safe_load(f)
    output = config['output']


def copy_files():
    '''Copy all files relevant for the clinic to the output folder'''

    df = pd.read_csv(f'{output}/samtools/minimum_error_rates.csv')
    for _, row in df.iterrows():
        read_id = row['read_id']
        ref = row['ref']

        consensus_in = f'{output}/consensus/medaka/{read_id}.{ref}/consensus.fasta'
        consensus_out = f'{output}/clinic/{read_id}.{ref}.fa'

        os.makedirs(f'{output}/clinic', exist_ok=True)

        # Change headers of consensus files and write to output folder
        with open(consensus_in, 'r') as in_fasta, open(consensus_out, 'w') as out_fa:
            for line in in_fasta:
                if line.startswith('>'):
                    out_fa.write(f'>{read_id}.{ref}\n')
                else:
                    out_fa.write(line)

    with open(f'{output}/copy_files_done.txt', 'w') as f:
        f.write('All files were copied to the output folder!')

copy_files()                                                                       