#!/usr/bin/env python3

import pandas as pd

df_files = pd.read_csv('samtools/minimum_error_rates.csv')    # File containing read_id and ref

for _, row in df_files.iterrows():
    
    # Get relevant vcf files 
    read_id = row['read_id']
    ref = row['ref']

    vcf_in_path = f'freebayes/{read_id}.{ref}.vcf'
    vcf_out_path = f'freebayes/{read_id}.{ref}_edit.vcf'

    # Read in vcf files and prepare for modifications
    with open(vcf_in_path, 'r') as f:
        vcf_in = f.readlines()
        header_line = next(line for line in vcf_in if line.startswith('#CHROM'))
        header = header_line.strip().split('\t')

        f.seek(0)
        
        df = pd.read_csv(f, comment='#', sep='\t', names=header)

    # Split the  column 'unknown' to extract genotype information
    if 'unknown' in df.columns and not df['unknown'].empty:
        split_df = df['unknown'].str.split(':', expand=True)                   
        split_df.columns = ['GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL']      
    else:
        print(f'Column "unknown" does not exist/is empty in {vcf_in_path}! Skipping...')
        continue

    df = pd.concat([df, split_df], axis=1)    # Concatenate the df with the split genotype columns
    df['GT'] = ' ' + df['GT']    # Avoid date conversion if vcf opened in Excel

    # Calculate the frequencies of mutations for variant calling
    split_ao = df['AO'].str.split(',', expand=True)
    df['RO'] = pd.to_numeric(df['RO'])

    for i in range(len(split_ao.columns)):
        ao_values = pd.to_numeric(split_ao[i])
        df[f'freq_{i+1:02}'] = (ao_values / (ao_values + df['RO'])).round(3)    # Calculate the frequencies (AO/(AO+RO)

    # Prepare data for output
    vcf_metadata = ''.join([line for line in vcf_in if line.startswith('##')])    # Lines containing metadata
    vcf_data = df.to_csv(sep='\t', index=False)    # Data lines
    output = vcf_metadata + vcf_data

    # Write frequencies to output vcf files
    with open(vcf_out_path, 'w') as f:
        f.write(output)
        with open('output/vcf_modifications_done.txt', 'w') as f:
            f.write('Modifications of vcf files done!')