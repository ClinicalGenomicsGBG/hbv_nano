#!/usr/bin/env python3

import pandas as pd

df = pd.read_csv("/medstore/Development/nanopore_HBV/daniel/test/error_rates.csv", header=None)    # import error rates

df.columns = ["read_id", "ref", "error_rate"]

df = df.loc[df.groupby("read_id")["error_rate"].idxmin()]   # get the sample with the lowest error rate for each read_id

print(df)