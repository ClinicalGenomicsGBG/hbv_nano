#!/usr/bin/env python3

import os
import gzip
import pandas as pd
'''
# QC of the negative control and sample(s)
def count_reads(fastq):    # read count of each fastq
    with gzip.open(fastq, 'rt') as f:
        num_lines = sum(1 for line in f)
    return num_lines // 4

fastq_files = [f for f in os.listdir('fastp') if f.endswith('_filt.fastq.gz')]    # get fastq files
samples = [f.replace('_filt.fastq.gz', '') for f in fastq_files]    # extract sample names from fastq files

reads = [count_reads(f'fastp/{f}') for f in fastq_files]    # count the number of reads in each fastq file

data = {'sample': samples, 'reads': reads}
df = pd.DataFrame(data)

reads_neg_ctrl = df.loc[df['sample'].str.contains('neg_ctrl'), 'reads'].values[0]
df['reads neg_ctrl/reads sample'] = round(reads_neg_ctrl / df['reads'], 3)

df['qc_pass_ctrl'] = df['reads neg_ctrl/reads sample'] <= 0.1
'''

'''
qc_check = df['qc_pass_ctrl'].sum() >= len(df) - 1   # check if all samples (except the negative control) has passed the QC

if qc_check == True:
    print("All samples have passed the QC")
else:
    print("ERROR: At least one sample has failed the QC")
'''
'''
# QC of coverage of the RT region
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
'''

##Under development##
df_error_rates = pd.read_csv('samtools/minimum_error_rates.csv')    # Use to select the samples
df_error_rates = df_error_rates[['read_id','ref']]

df_view_stats = pd.read_csv('samtools/view_stats.csv', header=None)
df_view_stats.columns = ['read_id', 'ref', 'mapped_reads']

# Read view_stats_rt.csv into a DataFrame
df_view_stats_rt = pd.read_csv('samtools/view_stats_rt.csv', header=None)
df_view_stats_rt.columns = ['read_id', 'ref', 'mapped_reads_rt']

# Merge the two DataFrames on read_id and ref
df_merged = pd.merge(df_error_rates, df_view_stats, on=['read_id', 'ref'])

df_merged = pd.merge(df_merged, df_view_stats_rt, on=['read_id', 'ref'])

reads_neg_ctrl = df_merged.loc[df_merged['read_id'].str.contains('neg_ctrl'), 'mapped_reads'].values[0]
reads_neg_ctrl_rt = df_merged.loc[df_merged['read_id'].str.contains('neg_ctrl'), 'mapped_reads_rt'].values[0]
#print(reads_neg_ctrl, reads_neg_ctrl_rt)
df_merged['mapped reads neg_ctrl/mapped reads sample'] = round(reads_neg_ctrl / df_merged['mapped_reads'], 3)
df_merged['mapped reads neg_ctrl_rt/mapped reads sample_rt'] = round(reads_neg_ctrl_rt / df_merged['mapped_reads_rt'], 3)

df_merged['qc_pass_ctrl'] = (df_merged['mapped reads neg_ctrl/mapped reads sample'] <= 0.1) & (df_merged['mapped_reads'] >= 100)   # Sample should have at least 10 times more reads than neg_ctrl and at least 100 mapped reads
df_merged['qc_pass_rt'] = (df_merged['mapped reads neg_ctrl_rt/mapped reads sample_rt'] <= 0.1) & (df_merged['mapped_reads_rt'] >= 100)    # Sample should have at least 10 times more reads than neg_ctrl and at least 100 reads mapped to the RT region 
df_merged['qc_pass'] = df_merged['qc_pass_ctrl'] & df_merged['qc_pass_rt']    # True if sample passes both QC checks
                                                                                                                                               
print(df_merged)

df_clinic = df_merged[['read_id','ref','mapped_reads','mapped_reads_rt','qc_pass']]

df_merged.to_csv('output/qc.csv', index=False)
df_clinic.to_csv('output/qc_clinic.csv', index=False)

##