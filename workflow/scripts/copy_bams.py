#!/usr/bin/env python3

import pandas as pd
import shutil
import os

def copy_bam_files():
    df = pd.read_csv('samtools/minimum_error_rates.csv')
    for index, row in df.iterrows():
        read_id = row['read_id']
        ref = row['ref']
        bam_file = f'samtools/{read_id}.{ref}.bam'
        if not os.path.exists('keep'):
            os.makedirs('keep')
        shutil.copy(bam_file, 'keep/')

copy_bam_files()