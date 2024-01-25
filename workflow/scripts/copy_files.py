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

        consensus_in = f'consensus/medaka/{read_id}.{ref}/consensus.fasta'
        consensus_out = f'output/{read_id}.{ref}_medaka.fa'
        
        #consensus_file = f'consensus/{read_id}.{ref}.fa'
        #medaka_dir = f'consensus/medaka/{read_id}.{ref}'
        #output_dir = f'output/medaka/{read_id}.{ref}'
        if not os.path.exists('output'):
            os.makedirs('output')

        with open(consensus_in, 'r') as in_fasta, open(consensus_out, 'w') as out_fa:
            for line in in_fasta:
                if line.startswith('>'):
                    out_fa.write(f'>{read_id}.{ref}\n')
                else:
                    out_fa.write(line)

        #shutil.copy(bam_file, 'output')
        #shutil.copy(consensus_file, 'output')
        #shutil.copytree(medaka_dir, output_dir)
        #shutil.copy(consensus_in, consensus_out)

    with open('output/copy_files_done.txt', 'w') as f:
        f.write('file copying done')

copy_files()                                                                       