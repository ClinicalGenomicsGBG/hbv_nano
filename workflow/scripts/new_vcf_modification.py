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

#output = get_output()    # Get the Snakemake output folder


def get_read_id_ref(file_path):
    '''Get the best matching reference for each read_id'''

    with open(file_path) as csv_file:
        reader = csv.DictReader(csv_file)
        return {row['read_id']: row['ref'] for row in reader}  #Check return!

#read_id_ref = get_read_id_ref(f'{output}/samtools/minimum_error_rates.csv')    # Read in file containing read_id and its reference
#print(f'{read_id_ref}, read_id_ref' )

def get_ref_genomes(ref_path):
    '''Get all the reference genomes (a->j))'''

    ref_genomes = {}

    for filepath in glob.glob(os.path.join(ref_path, 'ref_*.fa')):
        key = os.path.splitext(os.path.basename(filepath))[0]
        with open(filepath) as ref_path:
            ref_sequence = ''.join(line.strip() for line in ref_path if not line.startswith('>'))
        ref_genomes[key] = ref_sequence
    return ref_genomes


def read_vcf(reader):
    '''Read in the relevant information from the vcf file'''
    
    vcf = {}
    
    for record in reader:
        call = record.calls[0]
        pos = record.POS
        vcf[pos] = {
            'pos': pos,
            'ref': record.REF,
            'alt': [alt.value for alt in record.ALT] if record.ALT else [''],
            'qual': record.QUAL,
            'AO':  call.data.get('AO'),
            'RO': call.data.get('RO'),
        }
    return vcf  #TODO: döp om variabel till något annat

def split_vcf(vcf):
    '''Split the reference and and alternative sequnce(s) into sepparate positions in the vcf file.'''
    
    split_vcf = {}
    for pos, vcf in vcf.items():
        ref = vcf['ref']
        alts = vcf['alt']
        AO = vcf['AO']
        RO = vcf['RO']
        qual = vcf['qual']
        max_len = max([len(ref)] + [len(alt) for alt in alts])
        for i in range(max_len):
            key = pos + i
            split_row = {
                #'i': i,
                'pos': key,
                'ref': ref[i] if i < len(ref) else '',
                'qual': qual,
                'AO': AO,
                'RO': RO,
            }
            for idx, alt in enumerate(alts):
                split_row[f'alt_{idx + 1}'] = alt[i] if i < len(alt) else ''
                split_row[f'freq_{idx + 1}'] = round(AO[idx] / (sum(AO) + RO), 3)    # Calculate the frequency for each alt
            split_vcf[key] = split_row
    return split_vcf


def main():
    output = get_output()    # Get the Snakemake output folder
    ref_genomes = get_ref_genomes('reference_genomes')
    print(f'{ref_genomes}, ref_genomes')
    path = f'{output}/freebayes/KH20-2510.ref_d.vcf'
    reader = vcfpy.Reader.from_path(path)
    #vcf = vcf(reader)
    vcf = read_vcf(reader)
    vcf = split_vcf(vcf)
    '''
    for pos in range(130, 1161):    # The RT region
        if pos in vcf:
            print(vcf[pos])
    '''
if __name__ == '__main__':
    main()
