#!/usr/bin/env python3

import pandas as pd
import shutil
import os

def copy_files():
    df = pd.read_csv('samtools/minimum_error_rates.csv')
    for _, row in df.iterrows():
        read_id = row['read_id']
        ref = row['ref']
        #bam_file = f'samtools/{read_id}.{ref}.bam'
        consensus_file = f'consensus/{read_id}.{ref}.fa'
        medaka_dir = f'consensus/medaka/{read_id}.{ref}'
        output_dir = f'output/medaka/{read_id}.{ref}'
        if not os.path.exists('output'):
            os.makedirs('output')
        #shutil.copy(bam_file, 'output')
        shutil.copy(consensus_file, 'output')
        shutil.copytree(medaka_dir, output_dir)

    with open('output/copy_files_done.txt', 'w') as f:
        f.write('file copying done')

copy_files()                                                                       