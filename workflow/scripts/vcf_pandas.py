#!/usr/bin/env python3
#GENERALISING!!!
import pandas as pd

df_files = pd.read_csv('samtools/minimum_error_rates.csv')

for _, row in df_files.iterrows():
    read_id = row['read_id']
    ref = row['ref']

    vcf_in_path = f'freebayes/{read_id}.{ref}.vcf'
    vcf_out_path = f'freebayes/{read_id}.{ref}_edit.vcf'
    # Import vcf to df
    with open(vcf_in_path, 'r') as f:
        vcf_in = f.readlines()
        header_line = next(line for line in vcf_in if line.startswith('#CHROM'))
        header = header_line.strip().split('\t')

        f.seek(0)
        
        df = pd.read_csv(f, comment='#', sep='\t', names=header)

    if 'unknown' in df.columns and not df['unknown'].empty:
        split_df = df['unknown'].str.split(':', expand=True)                   # Split the 'unknown' column
        split_df.columns = ['GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL']    # and name the columns      
    else:
        print(f'Column "unknown" does not exist/is empty in {vcf_in_path}! Skipping...')
        continue

    df = pd.concat([df, split_df], axis=1)    # Concatenate the df with the split columns

    df['GT'] = ' ' + df['GT']    # Add space to avoid date conversion if vcf opened in Excel

    # Split the 'AO' column into multiple columns
    split_ao = df['AO'].str.split(',', expand=True)

    df['RO'] = pd.to_numeric(df['RO'])

    # Convert 'AO' column to list of numeric values and calculate AO/(AO+RO)
    for i in range(len(split_ao.columns)):
        ao_values = pd.to_numeric(split_ao[i])
        df[f'freq_{i+1:02}'] = (ao_values / (ao_values + df['RO'])).round(3)

    vcf_metadata = ''.join([line for line in vcf_in if line.startswith('##')])
    vcf_data = df.to_csv(sep='\t', index=False, header=False)
    output = vcf_metadata + vcf_data

    #print(vcf_data)
    # Write the DataFrame to the file
    with open(vcf_out_path, 'w') as f:
        f.write(output)

#print(df.filter(like='freq').iloc[1:11])