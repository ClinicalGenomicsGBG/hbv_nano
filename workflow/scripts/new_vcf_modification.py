#!/usr/bin/env python3

import pdb

import csv
import glob
import vcfpy
import os
import yaml


def get_output():
    '''Get output folder using snakemake or, if running script independently, directly from config file'''
    
    try:
        output = snakemake.params.output
    except NameError:
        with open('config/dev_config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        output = config['output']
    return output

output = get_output()    # Get the Snakemake output folder


def get_read_id_ref(file_path):
    '''Get the best matching reference for each read_id'''

    with open(file_path) as csv_file:
        reader = csv.DictReader(csv_file)
        return {row['read_id']: row['ref'] for row in reader}

read_id_ref = get_read_id_ref(f'{output}/samtools/minimum_error_rates.csv')    # Read in file containing read_id and its reference


def get_ref_genomes(ref_path):
    '''Get all the reference genomes (a->j))'''
    
    return {
        os.path.splitext(os.path.basename(filepath))[0]:
            ''.join(line.strip() for line in open(filepath) if not line.startswith('>'))
        for filepath in glob.glob(os.path.join(ref_path, 'ref_*.fa'))
    }

ref_genomes= get_ref_genomes('reference_genomes')


###dev
#reader = vcfpy.Reader.from_path('')
#record = next(reader)
#print(f'chrom: {record.CHROM}, pos: {record.POS}, ref: {record.REF}, alt: {record.ALT}')
###dev

def vcf(reader):
    vcf = {}
    for record in reader:
        call = record.calls[0]
        ref = record.REF
        alts = [alt.value for alt in record.ALT] if record.ALT else ['']
        pos = record.POS
        max_len = max([len(ref)] + [len(alt) for alt in alts])
        for i in range(max_len):
            key = pos + i
            split_row = {
                'pos': pos + i,
                'ref': ref[i] if i < len(ref) else '',
                'qual': record.QUAL,
                'AO': call.data.get('AO'),
                'RO': call.data.get('RO'),
            }
            for idx, alt in enumerate(alts):
                split_row[f'alt_{idx + 1}'] = alt[i] if i < len(alt) else ''
            vcf[key] = split_row
    return vcf


path = f'{output}/freebayes/KH20-2510.ref_d.vcf'
reader = vcfpy.Reader.from_path(path)
vcf = vcf(reader)

for pos in range(1700, 1900):
    if pos in vcf:
        print(vcf[pos])
