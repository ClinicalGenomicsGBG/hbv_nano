#!/usr/bin/env python3

import pdb

import csv
import yaml


def get_output():
    '''Get output folder using snakemake or, if running script independently, directly from config file'''
    
    try:
        output = snakemake.params.output
    except NameError:
        with open('config/config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        output = config['output']
    return output
breakpoint()
output = get_output()    # Get the Snakemake output folder
output = "/clinical/data/hbv_nano/results/250617-094351_250616_hbv_val_03_test"


def get_read_id_ref(file_path):
    '''Get the best matching reference for each read_id'''

    with open(file_path) as csv_file:
        reader = csv.DictReader(csv_file)
        return {row['read_id']: row['ref'] for row in reader}

read_id_ref = get_read_id_ref(f'{output}/samtools/minimum_error_rates.csv')    # Read in file containing read_id and its reference

