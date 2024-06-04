#!/usr/bin/env python3

import pandas as pd
import yaml

with open('config/config.yaml', 'r') as f:
    config = yaml.safe_load(f)
output = config['output']

# Calculate the minimum error rate for each read_id
#input_file = "samtools/error_rates.csv"
input_file = f'{output}/samtools/error_rates.csv'
#output_file = "output/samtools/minimum_error_rates.csv"
output_file = f'{output}/samtools/minimum_error_rates.csv'

df = pd.read_csv(input_file, header=None)    # import error rates
df.columns = ["read_id", "ref", "error_rate"]
df = df.loc[df.groupby("read_id")["error_rate"].idxmin()].reset_index(drop=True)   # get the sample with the lowest error rate for each read_id

df.to_csv(output_file, index=False)   # Write to csv