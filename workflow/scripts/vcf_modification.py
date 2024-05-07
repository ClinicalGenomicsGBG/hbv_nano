#!/usr/bin/env python3

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.Data import CodonTable


df_files = pd.read_csv('samtools/minimum_error_rates.csv')    # File containing read_id and ref

for _, row in df_files.iterrows():
    
    # Get relevant vcf files 
    read_id = row['read_id']
    ref = row['ref']

    vcf_in_path = f'freebayes/{read_id}.{ref}.vcf'
    vcf_out_path = f'freebayes/{read_id}.{ref}_edit.vcf'

    reference_genome = f'reference_genomes/{ref}.fa'
    
    with open(reference_genome, 'r') as f:
        lines = f.readlines()

    bases = ''.join(line.strip() for line in lines if not line.startswith('>'))
    ref_df = pd.DataFrame(list(bases), columns=['REF'])
    
    # Read in vcf files and prepare for modifications
    with open(vcf_in_path, 'r') as f:
        vcf_in = f.readlines()
        header_line = next(line for line in vcf_in if line.startswith('#CHROM'))
        header = header_line.strip().split('\t')

        f.seek(0)
        
        df = pd.read_csv(f, comment='#', sep='\t', names=header)

    # Split the  column 'unknown' to extract genotype information
    if 'unknown' in df.columns and not df['unknown'].empty:
        split_df = df['unknown'].str.split(':', expand=True)                   
        split_df.columns = ['GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL']      
    else:
        print(f'Column "unknown" does not exist/is empty in {vcf_in_path}! Skipping...')
        continue

    df = pd.concat([df, split_df], axis=1)    # Concatenate the df with the split genotype columns
    df['GT'] = ' ' + df['GT']    # Avoid date conversion if vcf opened in Excel
    
    # Split the column 'ALT' to separate columns ALT_01, ALT_02, etc.
    split_alt = df['ALT'].str.split(',', expand=True)
    
    for i in range(len(split_alt.columns)):
        split_alt.rename(columns={i: f'ALT_{i+1:02}'}, inplace=True)

    alt_index = df.columns.get_loc('ALT')
    df = df.drop(columns=['ALT'])    
    
    for i in reversed(range(len(split_alt.columns))):
        df.insert(alt_index, f'ALT_{i+1:02}', split_alt[f'ALT_{i+1:02}'])
    

    new_df = df[[col for col in df.columns if 'ALT' in col or col == 'QUAL']].copy()    # Create a new df with only the ALT columns (AND QUAL!!!)
     
    # Calculate the frequencies of mutations for variant calling
    split_ao = df['AO'].str.split(',', expand=True)
    df['RO'] = pd.to_numeric(df['RO'])  

    for i in range(len(split_ao.columns)):
        ao_values = pd.to_numeric(split_ao[i])
        new_df.loc[:,f'freq_{i+1:02}'] = (ao_values / (ao_values + df['RO'])).round(3)    # Calculate the frequencies (AO/(AO+RO)

                
    new_df.index = df['POS']    # Only get the positions where there are alternative alleles
    
    new_df.index.name = None
    
    out_df = ref_df    # Create output df starting with the reference df
    
    
    for col in [c for c in new_df.columns if 'ALT' in c or 'freq' in c]:    # Add the ALT columns to the output df
        out_df[col] = None                                                  # and fill with None
    
    out_df = ref_df.shift()    # Shift the reference dataframe to match the positions of the alternative alleles
    out_df['QUAL'] = None

    for idx, row in new_df.iterrows():    # Split the ALT columns into separate rows and put into the output df
        for col in [c for c in new_df.columns if 'ALT' in c]:
            if row[col] is not None:
                freq_col = col.replace('ALT', 'freq')
                for i, char in enumerate(row[col]):
                    out_df.loc[idx + i, col] = char
                    if freq_col in new_df.columns:    # Check if the 'freq' column exists
                        out_df.loc[idx + i, freq_col] = row[freq_col]
                        if 'QUAL' in new_df.columns:    # Check if the 'QUAL' column exists
                            out_df.loc[idx + i, 'QUAL'] = row['QUAL']
    
    # Extract the RT region
    start, end = 130, 1161    # Start and end of RT region
    positions = [(i, i+2) for i in range(start, end, 3)]    # The codons of the RT region
    
    new_df = pd.concat([out_df.loc[start:end] for start, end in positions])    # Extract the RT region from the output df

    alt_columns = [col for col in out_df.columns if 'ALT' in col]  # Get all 'ALT' columns

    for col in [c for c in new_df.columns if 'ALT' in c]:    # Fill empty values in the new df with the reference values
        new_df[col] = new_df[col].fillna(new_df['REF'])   
    
    for col in new_df.columns:    # Add empty 'aa_...' columns to the df
        new_df['aa_' + col] = None
    

    # Copy triplets of bases to the corresponding aa_columns
    for start, end in positions:
        for col in new_df.columns:
            if not col.startswith('aa_'):
                new_df.loc[start, 'aa_' + col] = new_df.loc[start:end, col].astype(str).str.cat()
            if 'QUAL' in new_df.columns:  # Check if the 'QUAL' column exists
                new_df.loc[start, 'aa_QUAL'] = new_df.loc[start:end, 'QUAL'].astype(str).str.cat(sep='; ')
            if col.startswith('freq'):  # Check if the column starts with 'freq'
                new_df.loc[start, 'aa_' + col] = new_df.loc[start:end, col].astype(str).str.cat(sep='; ')


    # If an ALT is equal to the REF remove it
    for col in [c for c in new_df.columns if c.startswith('aa_ALT')]:
        new_df.loc[new_df[col] == new_df['aa_REF'], col] = None

    if 'REF' in new_df.columns:
        new_df = new_df.drop(columns='REF')

    new_df = new_df.drop(columns=[col for col in new_df.columns if col.startswith('ALT_')])

 
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    amino_acid_codes = {
        'A': 'A - Alanine', 'R': 'R - Arginine', 'N': 'N - Asparagine', 'D': 'D - Aspartic acid', 'C': 'C - Cysteine',
        'E': 'E - Glutamic acid', 'Q': 'Q - Glutamine', 'G': 'G - Glycine', 'H': 'H - Histidine', 'I': 'I - Isoleucine',
        'L': 'L - Leucine', 'K': 'K - Lysine', 'M': 'M - Methionine', 'F': 'F - Phenylalanine', 'P': 'P - Proline',
        'S': 'S - Serine', 'T': 'T - Threonine', 'W': 'W - Tryptophan', 'Y': 'Y - Tyrosine', 'V': 'V - Valine',
        '*': 'Stop'
    }
    
    def translate_codon(codon):
        if codon is not None and codon in table.forward_table:
            return f"{amino_acid_codes[table.forward_table[codon]]} ({codon})"
        else:
            return None

    # Apply the translation function to the codon columns
    for col in [c for c in new_df.columns if 'aa_' in c and not 'freq' in c and c != 'aa_QUAL']:
        new_df[col] = new_df[col].apply(translate_codon)


    # Filter the dataframe
    alt_cols = [col for col in new_df.columns if 'aa_ALT' in col and not 'freq' in col]
    new_df = new_df[new_df['aa_REF'].notna() & new_df['aa_ALT_01'].notna() & (new_df[alt_cols].ne(new_df['aa_REF'], axis=0).any(axis=1))]    
    
    cols_to_drop = [col for col in new_df.columns if col.startswith('freq') or col == 'QUAL']
    new_df.drop(columns=cols_to_drop, inplace=True)
    
    new_df = new_df.replace('None', '__', regex = True)
    new_df.columns = new_df.columns.str.replace('aa_', '')    # Remove 'aa_' from the column names
    positions = [(367, 369), (634, 636), (646, 648), (667, 669), (670, 672),    # Known resistance positions
                (679,681), (733, 735), (739, 741), (835, 837), (877, 879)]
    
    new_df['Known resistance position'] = ''
    
    pd.options.display.float_format = '{:.0f}'.format
    new_df.insert(0, 'aa RT', ((new_df.index - 127) / 3).astype(int))    # Column with AA positions in RT region
    

    for start, end in positions:
        new_df.loc[start:end, 'Known resistance position'] = 'True'

    # Change triplets of idetical values to one entry to simplify readability
    def qual_simplify(value):
        parts = [part.strip() for part in value.split(';')]
        return parts [0] if len(set(parts)) == 1 else value
    
    new_df['QUAL'] = new_df['QUAL'].apply(qual_simplify)    # Change QUAL triplets to single value

    freq_cols = [col for col in new_df.columns if 'freq' in col]    # Extract all freq columns

    for col in freq_cols:   
        new_df[col] = new_df[col].apply(qual_simplify)    # Change freq triplets to single value

    
    ## Filter the output file based on the QUAL values
    new_df['qc_pass'] = ''

    # Function to check the QUAL values
    def qual_check(row):
        qual_val = row['QUAL'].split(';')
        qual_val = [float(x.strip()) for x in qual_val if x.strip() != '__']

        if not qual_val or all(x <= 1 for x in qual_val):
            return pd.Series(dtype='float64')

        if all(x <= 30 for x in qual_val):
            row['qc_pass'] = "False"
        
        if all(x > 30 for x in qual_val):
            row['qc_pass'] = "True"

        return row

    new_df = new_df.apply(qual_check, axis=1).dropna(how='all')    # Filter using the qual_check function
    
    # Move qc_pass column to the end of the df
    qc = new_df.pop('qc_pass')
    new_df = pd.concat([new_df, qc], axis=1)

    # Find the resistance mutations of the RT region 
    resistance_rows = pd.concat([new_df.loc[start:end] for start, end in positions])
    relevant_columns = ['aa RT', 'REF'] + alt_columns

    def amino_acid(row):
        aa_rt = str(int(row['aa RT']))
        for col in columns_to_modify:
            value = row[col]
            if pd.notna(value):
                row[col] = aa_rt + value[0]
        return row

    columns_to_modify = ['REF'] + alt_columns
    resistance_rows = resistance_rows.apply(amino_acid, axis=1)
    resistance_rows = resistance_rows[relevant_columns]
    print(resistance_rows, "resistance rows")
    

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

    for alt in alt_columns:
        # Check if the resistance mutations are present in the RT region
        for conditions, message in drugs.items():
            drug = conditions[-1]  # Get the drug name
            mutation_pairs = conditions[:-1]  # Get the mutation pairs
            present_mutations = [pair for pair in mutation_pairs if all(mutation in resistance_rows[alt].values for mutation in pair)]

            # Check for mutations and add write to output 
            if present_mutations: # if there are mutations for the alt

                present_mutations_str = ''.join(' and '.join(pair) for pair in present_mutations)
                category = message.split('; ')[0]
                data = pd.DataFrame([{'aa RT': present_mutations_str, 'REF': message}], columns=new_df.columns)
 
                if pd.notna(resistance_df.loc[resistance_df['Drug'] == drug, category]).any():
                    resistance_df.loc[resistance_df['Drug'] == drug, category] += ', ' + present_mutations_str
                else:
                    resistance_df.loc[resistance_df['Drug'] == drug, category] = present_mutations_str
                resistance_df.loc[resistance_df['Drug'] == drug, 'alt'] = alt

    col_names = resistance_df.columns.tolist()
    col_df = pd.DataFrame([col_names], columns=resistance_df.columns)
    empty_row_df = pd.DataFrame([dict.fromkeys(col_df.columns)], columns=col_df.columns)
    col_df = pd.concat([empty_row_df, col_df]).reset_index(drop=True)
    

    resistance_df = pd.concat([col_df, resistance_df]).reset_index(drop=True)

    resistance_df.columns = new_df.columns[:len(resistance_df.columns)]
    new_df = pd.concat([new_df, resistance_df], ignore_index=True)
    filename = f'output/variant_calling_{read_id}_{ref}.txt'
    new_df.to_csv(filename, sep='\t', index=True)

    with open('output/vcf_modifications_done.txt', 'w') as f:
        f.write('Modifications of vcf files done!')