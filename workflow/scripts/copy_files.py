#!/usr/bin/env python3

import pandas as pd
import shutil
import os
import yaml

with open('config/config.yaml', 'r') as f:
    config = yaml.safe_load(f)
output = config['output']


# Copy all files relevant for the clinic to the output folder
def copy_files():
    #df = pd.read_csv('samtools/minimum_error_rates.csv')
    df = pd.read_csv(f'{output}/samtools/minimum_error_rates.csv')
    for _, row in df.iterrows():
        read_id = row['read_id']
        ref = row['ref']

        consensus_in = f'consensus/medaka/{read_id}.{ref}/consensus.fasta'
        consensus_out = f'output/{read_id}.{ref}_medaka.fa'
        
        if not os.path.exists('output'):
            os.makedirs('output')

        # Change headers of consensus files and write to output folder
        with open(consensus_in, 'r') as in_fasta, open(consensus_out, 'w') as out_fa:
            for line in in_fasta:
                if line.startswith('>'):
                    out_fa.write(f'>{read_id}.{ref}\n')
                else:
                    out_fa.write(line)
        
        # Copy variant files to output folder
        #variant_in = f'freebayes/{read_id}.{ref}_edit.vcf'
        #variant_out = f'output/{read_id}.{ref}.vcf'
        #shutil.copy(variant_in, variant_out)

    #with open('output/copy_files_done.txt', 'w') as f:
    with open(f'{output}/copy_files_done.txt', 'w') as f:
        f.write('All files were copied to the output folder!')

copy_files()                                                                       