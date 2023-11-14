#!/usr/bin/env python3

"""
import pandas as pd

df = pd.read_csv("/medstore/Development/nanopore_HBV/daniel/test/error_rates.csv", header=None)    # import error rates

df.columns = ["read_id", "ref", "error_rate"]

df = df.loc[df.groupby("read_id")["error_rate"].idxmin()].reset_index(drop=True)   # get the sample with the lowest error rate for each read_id

#print(df)
df.to_csv("/medstore/Development/nanopore_HBV/daniel/test/error_rates_min.csv", index=False)   # export to csv
"""

import pandas as pd
import sys
import snakemake

#input_file = sys.argv[1]
#input_file = snakemake.input[0]
input_file = "/medstore/Development/nanopore_HBV/daniel/hbv_nano/samtools/error_rates.csv"

#output_file = sys.argv[2]
#output_file = snakemake.output[0]
output_file = "/medstore/Development/nanopore_HBV/daniel/hbv_nano/samtools/minimum_error_rates.csv"

#with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
#    pass
df = pd.read_csv(input_file, header=None)    # import error rates

df.columns = ["read_id", "ref", "error_rate"]

df = df.loc[df.groupby("read_id")["error_rate"].idxmin()].reset_index(drop=True)   # get the sample with the lowest error rate for each read_id

df.to_csv(output_file, index=False)   # export to csv