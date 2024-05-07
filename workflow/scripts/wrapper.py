#!/usr/bin/env python3

import pandas as pd
import glob
import shutil

sample_sheet = pd.read_csv('data/samplesheet.csv')

#print(sample_sheet)

def concatenate(sample_sheet, rootdir, outputdir):
    print("Concatenating files")    
    for _, row in sample_sheet.iterrows():
        alias = row['alias']
        bc_files = glob.glob(f"{rootdir}/{alias}/*{alias}*.fastq.gz")
        with open(f"{outputdir}/{alias}.fastq.gz",'wb') as wfd:
            for f in bc_files:
                with open(f,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)


concatenate(sample_sheet, "data/test_data", "data/out_data")