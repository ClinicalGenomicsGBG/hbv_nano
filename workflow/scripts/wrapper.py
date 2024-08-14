#!/usr/bin/env python3

import pandas as pd
import glob
import shutil
import os
import yaml

with open('config/wrapper_config.yaml', 'r') as f:
    config = yaml.safe_load(f)

# Get paths from config file
inpath = config['input_dir']    # Input path containing barcode folders with fastq files
outpath = config['output_dir']    # Output path for concatenated fastq files
sample_sheet = pd.read_csv(config['sample_sheet'])    # Sample sheet following Nanopore convention

# Create output directory if it does not exist (but is defined in the config file)
if not os.path.exists(outpath):
    os.makedirs(outpath)

# Check if there is a negative control in the sample sheet
if not sample_sheet['type'].str.contains('neg', na=False).any():
    raise ValueError("No negative control could be found, please check input data. Aborting script.")

# Function to concatenate fastq files for each barcode
def concatenate(sample_sheet, inpath, outpath):
    
    for _, row in sample_sheet.iterrows():
        barcode = row['barcode']
        alias = row['alias']
        sample_type = str(row['type']) if not pd.isna(row['type']) else ''
        fastq_files = glob.glob(f"{inpath}/{barcode}/*{barcode}*.fastq.gz")
        
        if not fastq_files:
            print(f"No files found for alias {barcode}")
            continue
            
        if 'neg' in sample_type:
            alias = f'{alias}_neg_ctrl'    # Change alias for negative control sample

        print(f'Concatenating all files for {barcode}.fastq.gz to {alias}.fastq.gz')
        
        with open(f"{outpath}/{alias}.fastq.gz",'wb') as wfd:
            for f in fastq_files:
                with open(f,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)

# Concatenate fastq files
concatenate(sample_sheet, inpath, outpath)
