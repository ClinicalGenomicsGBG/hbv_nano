#!/usr/bin/env python3

import pandas as pd

input_file = "samtools/error_rates.csv"
output_file = "samtools/minimum_error_rates.csv"
output_file_2 = "samtools/minimum_error_rates_2.csv"

df = pd.read_csv(input_file, header=None)    # import error rates
df.columns = ["read_id", "ref", "error_rate"]
df = df.loc[df.groupby("read_id")["error_rate"].idxmin()].reset_index(drop=True)   # get the sample with the lowest error rate for each read_id

df.to_csv(output_file, index=False)   # export to csv

new_df = df["read_id"] + "." + df["ref"] + ".bam"   # df with all the relevant bam files.
new_df.to_csv(output_file_2, index=False, header=False)