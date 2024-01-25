#!/usr/bin/env python3

import os
import gzip
import pandas as pd

def count_reads(fastq):
    with gzip.open(fastq, 'rt') as f:
        num_lines = sum(1 for line in f)
    return num_lines // 4

fastq_files = [f for f in os.listdir('fastp') if f.endswith('_filt.fastq.gz')]

samples = [f.replace('_filt.fastq.gz', '') for f in fastq_files]    # get filenames
reads = [count_reads(f'fastp/{f}') for f in fastq_files]

data = {'sample': samples, 'reads': reads}

df = pd.DataFrame(data)

reads_neg_ctrl = df.loc[df['sample'].str.contains('neg_ctrl'), 'reads'].values[0]
df['reads neg_ctrl/reads sample'] = round(reads_neg_ctrl / df['reads'], 3)

df['qc_pass_ctrl'] = df['reads neg_ctrl/reads sample'] <= 0.1    # check if the reads for the samples are less or equal 10% of the negative control reads

qc_check = df['qc_pass_ctrl'].sum() >= len(df) - 1   # check if all samples (except the negative control) has passed the QC

if qc_check == True:
    print("All samples have passed the QC")
else:
    print("ERROR: At least one sample has failed the QC")

consensus_in = [f for f in os.listdir('output') if f.endswith('_medaka.fa')]

for filename in consensus_in:
    with open(os.path.join('output', filename), 'r') as in_fasta:
        lines = in_fasta.readlines()
        sequence = ''.join(line.strip() for line in lines[1:])
        rt_region = sequence[129:1161]
        count = round(100 - (rt_region.count('N') / len(rt_region)) * 100, 2)
        sample_name = filename.split('.')[0]
        count_n = round(100 - (sequence.count('N')/len(sequence)) * 100, 2)
        df.loc[df['sample'] == sample_name, 'RT region coverage [%]'] = count
        df.loc[df['sample'] == sample_name, 'covered bases [%]'] = count_n

df['qc_pass_RT'] = df['RT region coverage [%]'] >= 90

df.insert(1, 'covered bases [%]', df.pop('covered bases [%]'))

df = df.sort_values('sample')
df.to_csv('output/qc.csv', index=False)