#!/usr/bin/env python3

import pandas as pd

# Get the best matching samples
df_error_rates = pd.read_csv('samtools/minimum_error_rates.csv')
df_error_rates = df_error_rates[['read_id','ref']]

# Get the number of mapped reads for each sample
df_view_stats = pd.read_csv('samtools/view_stats.csv', header=None)
df_view_stats.columns = ['read_id', 'ref', 'mapped_reads']

# Get the number of mapped reads in the rt region for each sample
df_view_stats_rt = pd.read_csv('samtools/view_stats_rt.csv', header=None)
df_view_stats_rt.columns = ['read_id', 'ref', 'mapped_reads_rt']

# Get the number of mapped reads for the best matching samples
df_merged = pd.merge(df_error_rates, df_view_stats, on=['read_id', 'ref'])
df_merged = pd.merge(df_merged, df_view_stats_rt, on=['read_id', 'ref'])

# Get the number of mapped reads for the negative controls of the samples and the rt regions
reads_neg_ctrl = df_merged.loc[df_merged['read_id'].str.contains('neg_ctrl'), 'mapped_reads'].values[0]
reads_neg_ctrl_rt = df_merged.loc[df_merged['read_id'].str.contains('neg_ctrl'), 'mapped_reads_rt'].values[0]

# Calculate the ratio of mapped reads of the negative control to the sample and the rt region
df_merged['mapped reads neg_ctrl/mapped reads sample'] = round(reads_neg_ctrl / df_merged['mapped_reads'], 3)
df_merged['mapped reads neg_ctrl_rt/mapped reads sample_rt'] = round(reads_neg_ctrl_rt / df_merged['mapped_reads_rt'], 3)

# Perform QC checks
df_merged['qc_pass_ctrl'] = (df_merged['mapped reads neg_ctrl/mapped reads sample'] <= 0.1) & (df_merged['mapped_reads'] >= 100)   # Sample should have at least 10 times more reads than neg_ctrl and at least 100 mapped reads
df_merged['qc_pass_rt'] = (df_merged['mapped reads neg_ctrl_rt/mapped reads sample_rt'] <= 0.1) & (df_merged['mapped_reads_rt'] >= 100)    # Sample should have at least 10 times more reads than neg_ctrl and at least 100 reads mapped to the RT region 
df_merged['qc_pass'] = df_merged['qc_pass_ctrl'] & df_merged['qc_pass_rt']    # True if sample passes both QC checks
                                                                                                                                               
print(df_merged)

# Output the results
# All results
df_merged.to_csv('output/qc.csv', index=False)

# Results relevant to the clinic
df_clinic = df_merged[['read_id','ref','mapped_reads','mapped_reads_rt','qc_pass']]
df_clinic.to_csv('output/qc_clinic.csv', index=False)