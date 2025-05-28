#!/usr/bin/env python3

import pandas as pd
import os
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.Data import CodonTable
import yaml
import csv

def get_output():
    '''Get output folder using snakemake or, if running script independently, directly from config file'''
    
    try:
        output = snakemake.params.output
    except NameError:
        with open('config/config.yaml', 'r') as f:
            config = yaml.safe_load(f)
        output = config['output']
    return output

output = get_output()    # Get the output folder path

def translate_codon(codon):
    '''Translate codons to amino acids'''

    if codon is not None and codon in table.forward_table:
        return f"{amino_acid_codes[table.forward_table[codon]]} ({codon})"
    else:
        return None

def get_samples(file_path):
    '''Get the best matchings samples'''

    sample_dict = {}
    with open(file_path) as csv_file:
        reader = csv.DictReader(csv_file)
        for row in reader:
            key = row['read_id']
            sample_dict[key] = row['ref']
    return sample_dict

min_err_rates = get_samples(f'{output}/samtools/minimum_error_rates.csv')
samples = get_samples(f'{output}/samtools/minimum_error_rates.csv')    # Read in file containing read_id and ref


def get_references(reference_path):
    '''Read all reference genome files in a directory and return a dictionary with ref values as keys and sequences as values.'''
    reference_dict = {}
    valid_files = {'ref_a.fa', 'ref_b.fa', 'ref_c.fa', 'ref_d.fa', 'ref_e.fa', 'ref_f.fa', 'ref_g.fa', 'ref_h.fa', 'ref_i.fa', 'ref_j.fa'}
    for reference in os.listdir(reference_path):
        if reference in valid_files:
            reference_genome = os.path.join(reference_path, reference)
            base_name, _ = os.path.splitext(reference)
            with open(reference_genome, 'r') as f:
                next(f)  # Skip the header line
                sequence = ''.join(line.strip() for line in f)
                reference_dict[base_name] = sequence
    return reference_dict

reference_genomes = get_references('reference_genomes')    # Read in reference genomes
#####



for read_id, ref in samples.items():    

    vcf_in_path = f'{output}/freebayes/{read_id}.{ref}.vcf'

    reference_genome = f'reference_genomes/{ref}.fa'
    
    # Extract sequences from reference genomes
    with open(reference_genome, 'r') as f:
        lines = f.readlines()

    bases = ''.join(line.strip() for line in lines if not line.startswith('>'))
    ref_df = pd.DataFrame(list(bases), columns=['REF'])
    
    
    # Read in vcf file and prepare for modifications
    with open(vcf_in_path, 'r') as f:
        vcf_in = f.readlines()
        header_line = next(line for line in vcf_in if line.startswith('#CHROM'))
        header = header_line.strip().split('\t')
        
        f.seek(0)    # Go back to the beginning of the file
        
        vc_df = pd.read_csv(f, comment='#', sep='\t', names=header)    # Read in the vcf file using header as column names
        
    # Split the  column 'unknown' to extract genotype information
    if 'unknown' in vc_df.columns and not vc_df['unknown'].empty:
        split_df = vc_df['unknown'].str.split(':', expand=True)                   
        split_df.columns = ['GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL']      
    else:
        print(f'Column "unknown" does not exist/is empty in {vcf_in_path}! Skipping...')
        continue

    vc_df = pd.concat([vc_df, split_df], axis=1)    # Concatenate the df with the split genotype columns
    vc_df['GT'] = ' ' + vc_df['GT']    # Avoid date conversion if vcf opened in Excel
    
    # Split the column 'ALT' to separate columns ALT_01, ALT_02, etc.
    split_alt = vc_df['ALT'].str.split(',', expand=True)
    
    for i in range(len(split_alt.columns)):
        split_alt.rename(columns={i: f'ALT_{i+1:02}'}, inplace=True)

    alt_index = vc_df.columns.get_loc('ALT')
    vc_df = vc_df.drop(columns=['ALT'])    
    
    for i in reversed(range(len(split_alt.columns))):
        vc_df.insert(alt_index, f'ALT_{i+1:02}', split_alt[f'ALT_{i+1:02}'])
    
    aa_data_df = vc_df[[col for col in vc_df.columns if 'ALT' in col or col == 'QUAL']].copy()    # df containing the ALT columns and the QUAL column

    # Calculate the frequencies of mutations for variant calling
    split_ao = vc_df['AO'].str.split(',', expand=True)
    split_ao = split_ao.apply(pd.to_numeric).fillna(0)
    ao_tot = split_ao.sum(axis=1)
    
    vc_df['RO'] = pd.to_numeric(vc_df['RO'])

    for i in range(len(split_ao.columns)):
        ao_values = split_ao[i]
        aa_data_df.loc[:,f'freq_{i+1:02}'] = (ao_values / (ao_tot + vc_df['RO'])).round(3)
                
    aa_data_df.index = vc_df['POS']    # Only get the positions where there are alternative alleles
    aa_data_df.index.name = None
    
    out_df = ref_df    # Create output df starting with the reference df
    
    
    for col in [c for c in aa_data_df.columns if 'ALT' in c or 'freq' in c]:    # Add the ALT columns to the output df
        out_df[col] = None                                                      # and fill with None
    
    out_df = ref_df.shift()    # Shift the reference dataframe to match the positions of the alternative alleles
    out_df['QUAL'] = None


    for idx, row in aa_data_df.iterrows():    # Split the ALT columns into separate rows and put into the output df
        for col in [c for c in aa_data_df.columns if 'ALT' in c]:
            if row[col] is not None:
                
                freq_col = col.replace('ALT', 'freq')    # Add freq columns. E.g. ALT_01 -> freq_01...
                
                for i, char in enumerate(row[col]):
                    out_df.loc[idx + i, col] = char    # Takes the alternative sequences from the columns, splits them and puts them into the output df row by row
                    #print(f'idx: {idx}, i: {i}, col: {col}, char: {char}') #test
                    if freq_col in aa_data_df.columns:    # Check if the 'freq' column exists
                        out_df.loc[idx + i, freq_col] = row[freq_col]    # Add the frequency values to out_df
                    if 'QUAL' in aa_data_df.columns:    # Check if the 'QUAL' column exists
                        out_df.loc[idx + i, 'QUAL'] = row['QUAL']    # Add the QUAL values to out_df
    
    # Extract the RT region
    start, end = 130, 1161    # Start and end of RT region
    positions = [(i, i+2) for i in range(start, end, 3)]    # The codons of the RT region
    
    aa_data_df = pd.concat([out_df.loc[start:end] for start, end in positions])    # Extract the RT region from the output df

    alt_columns = [col for col in out_df.columns if 'ALT' in col]  # Get all 'ALT' columns

    for col in [c for c in aa_data_df.columns if 'ALT' in c]:    # Fill empty values in the new df with the reference values
        aa_data_df[col] = aa_data_df[col].fillna(aa_data_df['REF'])   
    
    for col in aa_data_df.columns:    # Add empty 'aa_...' columns to the df
        aa_data_df['aa_' + col] = None
    

    # Copy triplets of bases to the corresponding aa_columns
    for start, end in positions:
        for col in aa_data_df.columns:
            if not col.startswith('aa_'):
                aa_data_df.loc[start, 'aa_' + col] = aa_data_df.loc[start:end, col].astype(str).str.cat()
            if 'QUAL' in aa_data_df.columns:  # Check if the 'QUAL' column exists
                aa_data_df.loc[start, 'aa_QUAL'] = aa_data_df.loc[start:end, 'QUAL'].astype(str).str.cat(sep='; ')
            if col.startswith('freq'):  # Check if the column starts with 'freq'
                aa_data_df.loc[start, 'aa_' + col] = aa_data_df.loc[start:end, col].astype(str).str.cat(sep='; ')


    # If an ALT is equal to the REF remove it
    for col in [c for c in aa_data_df.columns if c.startswith('aa_ALT')]:
        aa_data_df.loc[aa_data_df[col] == aa_data_df['aa_REF'], col] = None

    if 'REF' in aa_data_df.columns:
        aa_data_df = aa_data_df.drop(columns='REF')

    aa_data_df = aa_data_df.drop(columns=[col for col in aa_data_df.columns if col.startswith('ALT_')])

    table = CodonTable.unambiguous_dna_by_name["Standard"]
    amino_acid_codes = {
        'A': 'A - Alanine', 'R': 'R - Arginine', 'N': 'N - Asparagine', 'D': 'D - Aspartic acid', 'C': 'C - Cysteine',
        'E': 'E - Glutamic acid', 'Q': 'Q - Glutamine', 'G': 'G - Glycine', 'H': 'H - Histidine', 'I': 'I - Isoleucine',
        'L': 'L - Leucine', 'K': 'K - Lysine', 'M': 'M - Methionine', 'F': 'F - Phenylalanine', 'P': 'P - Proline',
        'S': 'S - Serine', 'T': 'T - Threonine', 'W': 'W - Tryptophan', 'Y': 'Y - Tyrosine', 'V': 'V - Valine',
        '*': 'Stop'
    }

    # Apply the translation function to the codon columns
    for col in [c for c in aa_data_df.columns if 'aa_' in c and not 'freq' in c and c != 'aa_QUAL']:
        aa_data_df[col] = aa_data_df[col].apply(translate_codon)


    # Filter the dataframe
    alt_cols = [col for col in aa_data_df.columns if 'aa_ALT' in col and not 'freq' in col]
    aa_data_df = aa_data_df[aa_data_df['aa_REF'].notna() & aa_data_df['aa_ALT_01'].notna() & (aa_data_df[alt_cols].ne(aa_data_df['aa_REF'], axis=0).any(axis=1))]    
    
    cols_to_drop = [col for col in aa_data_df.columns if col.startswith('freq') or col == 'QUAL']
    aa_data_df.drop(columns=cols_to_drop, inplace=True)
    
    aa_data_df = aa_data_df.replace('None', '__', regex = True)
    aa_data_df.columns = aa_data_df.columns.str.replace('aa_', '')    # Remove 'aa_' from the column names
    positions = [(367, 369), (634, 636), (646, 648), (667, 669), (670, 672),    # Known resistance positions
                (679,681), (733, 735), (739, 741), (835, 837), (877, 879)]
    
    aa_data_df['Known resistance position'] = ''
    
    pd.options.display.float_format = '{:.0f}'.format
    aa_data_df.insert(0, 'aa RT', ((aa_data_df.index - 127) / 3).astype(int))    # Column with AA positions in RT region
    

    for start, end in positions:
        aa_data_df.loc[start:end, 'Known resistance position'] = 'True'

    def qual_simplify(value):
        '''Change triplets of idetical values to one entry to simplify readability'''

        parts = [part.strip() for part in value.split(';')]
        return parts [0] if len(set(parts)) == 1 else value
    
    aa_data_df['QUAL'] = aa_data_df['QUAL'].apply(qual_simplify)    # Change QUAL triplets to single value

    freq_cols = [col for col in aa_data_df.columns if 'freq' in col]    # Extract all freq columns

    for col in freq_cols:   
        aa_data_df[col] = aa_data_df[col].apply(qual_simplify)    # Change freq triplets to single value

    
    ## Filter the output file based on the QUAL values
    aa_data_df['qc_pass'] = ''

    def qual_check(row):
        '''Check the QUAL values'''

        qual_val = row['QUAL'].split(';')
        qual_val = [float(x.strip()) for x in qual_val if x.strip() != '__']

        # If all qual values are either non-existent or <= 1 and assign NA
        if not qual_val or all(x <= 1 for x in qual_val):
            row['qc_pass'] = "NA"
            return row
        
        # Check if any qual value <= 30 and assign True, False
        if any(x <= 30 for x in qual_val):
            row['qc_pass'] = "False"
        else:
            row['qc_pass'] = "True"
        return row

    # Filter based on the QUAL values
    aa_data_df = aa_data_df.apply(qual_check, axis=1).dropna(how='all')    # Filter using the qual_check function
    aa_data_df = aa_data_df[aa_data_df['qc_pass'] != "NA"]    # Drop rows where all qual values are either non-existent or <= 1
    
    # Move qc_pass column to the end of the df
    qc = aa_data_df.pop('qc_pass')
    aa_data_df = pd.concat([aa_data_df, qc], axis=1)

    # Find the resistance mutations of the RT region 
    resistance_rows = pd.concat([aa_data_df.loc[start:end] for start, end in positions])
    relevant_columns = ['aa RT', 'REF'] + alt_columns

    def amino_acid(row):
        '''Modify the amino acid columns to include the position in the RT region'''
        
        aa_pos_rt = str(int(row['aa RT']))    # Get the AA position in the RT region
        
        for col in aa_columns:    # Modify the REF column and the ALT column(s)
            aa_name = row[col]
            if pd.notna(aa_name):
                row[col] = aa_pos_rt + aa_name[0]    # Add the RT position to the AA name
        return row

    aa_columns = ['REF'] + alt_columns
    resistance_rows = resistance_rows.apply(amino_acid, axis=1)
    resistance_rows = resistance_rows[relevant_columns]
    

    # Amino acid changes and the drug resitance they cause
    drugs = {
        (("173L",), "Lamivudine"): "Compensatory mutation; Lamivudine",
        (("180C",), ("180M",), "Lamivudine"): "Limited susceptibility; Lamivudine",
        (("204I",), ("204S",), ("204V",), "Lamivudine"): "Resistant; Lamivudine",
        (("181T",), ("181V",), "Lamivudine"): "Resistant; Lamivudine",
        (("80V"), ("80I",), "Lamivudine"): "Compensatory mutation; Lamivudine",

        (("181T",), ("181V",), "Adefovir"): "Resistant; Adefovir, Hepsera",      
        (("236T",), "Adefovir"): "Resistant; Adefovir, Hepsera",

        (("169T","204V"), ("184A","204V"), ("184G","204V"), ("184I","204V"), ("184S","204V"), ("202G","204V"), ("202I","204V"),("250V","204V"), "Entecavir"): "Resistant; Entecavir",
        (("169T","204I"), ("184A","204I"), ("184G","204I"), ("184I","204I"), ("184S","204I"), ("202G","204I"), ("202I","204I"), ("250V","204I"), "Entecavir"): "Resistant; Entecavir",

        (("204V",), ("204I",), "Entecavir"): "Partly resistant; Entecavir",
        (("180C",), ("180M",), "Entecavir"): "Compensatory mutation; Entecavir",
        
        (("169T",), ("184A",), ("184G",), ("184I",), ("184S",), "Entecavir"): "Compensatory mutation; Entecavir",
        (("202G",), ("202I",), "Entecavir"): "Compensatory mutation; Entecavir",
        (("250V",), "Entecavir"): "Compensatory mutation; Entecavir",

        (("204I",), ("204V",), "Telbivudine"): "Resistant; Telbivudine",
        (("80I",), ("80V",), "Telbivudin"): "Compensatory mutation; Telbivudine",
        (("181T",), ("181V",), "Telbivudine"): "Resistant; Telbivudine",
    }
    
    resistance_df = pd.DataFrame(columns=['Drug', 'Compensatory mutation', 'Limited susceptibility', 'Partly resistant', 'Resistant', 'alt'])
    drug_list = ["Entecavir", "Lamivudine", "Telbivudine", "Adefovir", "Tenofovir"]
    resistance_df['Drug'] = drug_list
    
    resistance_dfs = {} # Store AA changes for each ALT
    all_dfs = [] # Store all dfs for each ALT
    

    for alt in alt_columns:
        resistance_dfs[alt] = resistance_df.copy()
        
        # Check if the resistance mutations are present in the RT region
        for conditions, message in drugs.items():
            drug = conditions[-1]  # Get the drug name
            mutation_pairs = conditions[:-1]  # Get the mutation pairs
            present_mutations = [pair for pair in mutation_pairs if all(mutation in resistance_rows[alt].values for mutation in pair)]

            # Check for mutations and add write to output 
            if present_mutations: # if there are mutations for the alt

                present_mutations_str = ''.join(' and '.join(pair) for pair in present_mutations)
                category = message.split('; ')[0]
                data = pd.DataFrame([{'aa RT': present_mutations_str, 'REF': message}], columns=aa_data_df.columns)

                if pd.notna(resistance_dfs[alt].loc[resistance_dfs[alt]['Drug'] == drug, category]).any():
                    resistance_dfs[alt].loc[resistance_dfs[alt]['Drug'] == drug, category] += f', {present_mutations_str} ({alt})'
                    resistance_dfs[alt].loc[resistance_dfs[alt]['Drug'] == drug, 'alt'] = alt
                else:
                    resistance_dfs[alt].loc[resistance_dfs[alt]['Drug'] == drug, category] = f'{present_mutations_str} ({alt})'
                    resistance_dfs[alt].loc[resistance_dfs[alt]['Drug'] == drug, 'alt'] = alt

                resistance_df.loc[resistance_df['Drug'] == drug]

        all_dfs.append(resistance_dfs[alt])

    ## Concatenate the drug resistance for each alt and format for output
    output_resistance_df = pd.concat([pd.concat([vc_df, pd.DataFrame([{}])], ignore_index=True) for vc_df in all_dfs], ignore_index=True)
    header_resistance_df = pd.DataFrame(columns = resistance_df.columns)
    header_resistance_df.loc[0] = resistance_df.columns
    resistance_df = pd.concat([header_resistance_df, output_resistance_df], ignore_index=True)
    resistance_df.columns = aa_data_df.columns[:len(resistance_df.columns)]
    empty_row_df = pd.DataFrame([{}], columns=aa_data_df.columns)
    output_df = pd.concat([aa_data_df, empty_row_df, resistance_df], ignore_index=True)    # Final output df

    # Finalise and write output
    filename = f'{output}/variant_calling_{read_id}_{ref}.txt'
    
    with open(filename, 'w') as f:
        output_df.to_csv(filename, sep='\t', index=True)

    with open(f'{output}/vcf_modifications_done.txt', 'w') as f:
        f.write('Modifications of vcf files done!')
    