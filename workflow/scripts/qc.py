#!/usr/bin/env python3

import pandas as pd
import yaml
import os


def get_output():
    '''Get output folder using snakemake or, if running script independently, directly from config file'''
    
    try:
        output = snakemake.params.output
    except NameError:
        with open('config/config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        output = config['output']
    return output

def data_handling():
    '''Get mapped reads best matching samples and negative controls and calculate ratios'''

    output = get_output()

    # Get the best matching samples
    df_min_error_rates = pd.read_csv(f'{output}/samtools/minimum_error_rates.csv')
    df_error_rates = df_min_error_rates[['read_id','ref']]

    # Get the number of mapped reads for each sample
    df_view_stats = pd.read_csv(f'{output}/samtools/view_stats.csv', header=None)
    df_view_stats.columns = ['read_id', 'ref', 'mapped_reads']

    # Get the number of mapped reads in the rt region for each sample
    df_view_stats_rt = pd.read_csv(f'{output}/samtools/view_stats_rt.csv', header=None)
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

    return df_merged

def qc(df_merged):
    '''Perform QC checks'''

    # Perform QC checks
    df_merged['qc_pass_ctrl'] = (df_merged['mapped reads neg_ctrl/mapped reads sample'] <= 0.1) & (df_merged['mapped_reads'] >= 100)   # Sample should have at least 10 times more reads than neg_ctrl and at least 100 mapped reads
    df_merged['qc_pass_rt'] = (df_merged['mapped reads neg_ctrl_rt/mapped reads sample_rt'] <= 0.1) & (df_merged['mapped_reads_rt'] >= 100)    # Sample should have at least 10 times more reads than neg_ctrl and at least 100 reads mapped to the RT region 
    df_merged['qc_pass'] = df_merged['qc_pass_ctrl'] & df_merged['qc_pass_rt']    # True if sample passes both QC checks

    # Perform QC check of the negative control
    mapped_reads_tot = df_merged['mapped_reads'].sum()
    reads_neg_ctrl = df_merged.loc[df_merged['read_id'].str.contains('neg_ctrl'), 'mapped_reads'].values[0]
    neg_ctrl_ratio = round(reads_neg_ctrl / mapped_reads_tot, 3)

    if neg_ctrl_ratio > 0.1:
        df_merged.loc[df_merged['read_id'].str.contains('neg_ctrl'), 'qc_pass'] = False
        print(f'Warning! The negative control contains {neg_ctrl_ratio * 100}% of the reads')
    else:
        df_merged.loc[df_merged['read_id'].str.contains('neg_ctrl'), 'qc_pass'] = True

    return df_merged

def output_results(df_merged, output):
    '''Output the results'''

    # Output full qc file
    df_merged.to_csv(f'{output}/qc_full.csv', index=False)

    # Output qc file with results relevant to the clinic
    df_clinic = df_merged[['read_id','ref','mapped_reads','mapped_reads_rt','qc_pass']]
    df_clinic.to_csv(f'{output}/qc.csv', index=False)

def main():
    output = get_output()
    df_merged = data_handling()
    df_merged = qc(df_merged)
    output_results(df_merged, output)

if __name__ == '__main__':
    main()
