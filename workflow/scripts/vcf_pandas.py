#!/usr/bin/env python3

import pandas as pd

# Read the header lines
with open('freebayes/bc17.ref_d.vcf', 'r') as f:
    lines = f.readlines()
    header_line = next(line for line in lines if line.startswith('#CHROM'))
    header = header_line.strip().split('\t')

# Read the data lines into a DataFrame
df = pd.read_csv('freebayes/bc17.ref_d.vcf', comment='#', sep='\t', names=header)

# Split the 'unknown' column into multiple columns
split_df = df['unknown'].str.split(':', expand=True)
split_df.columns = ['GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL']

# Concatenate the original DataFrame with the new DataFrame
df = pd.concat([df, split_df], axis=1)

df['GT'] = ' ' + df['GT']

#print(df)
# Split the 'AO' column into multiple columns
split_ad = df['AO'].str.split(',', expand=True)

# Convert 'RO' column to numeric
df['RO'] = pd.to_numeric(df['RO'])

# Convert 'AO' column to list of numeric values and calculate AO/(AO+RO)
for i in range(len(split_ad.columns)):
    ao_values = pd.to_numeric(split_ad[i])
    print(ao_values)
    df[f'freq_{i+1:02}'] = (ao_values / (ao_values + df['RO'])).round(3)

metadata = ''.join([line for line in lines if line.startswith('##')])
df_str = df.to_csv(sep='\t', index=False, header=True)
output = metadata + df_str

# Print the 'freq' columns
print(df.filter(like='freq'))

# Write the DataFrame to the file
with open('freebayes/bc17.ref_d_edit_pandassssss.vcf', 'w') as f:
    f.write(output)

#print(df)