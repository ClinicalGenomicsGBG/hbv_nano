# Reference sequences (a-i)

ref_folder = "/Users/xschmd/Desktop/referensgenom"    #location of reference sequences

#print(glob.glob(f"{ref_folder}/*.fa"))

ref_a = f"{ref_folder}/ref_a.fa"    #location of reference sequence
ref_b = f"{ref_folder}/ref_b.fa"
ref_c = f"{ref_folder}/ref_c.fa"
ref_d = f"{ref_folder}/ref_d.fa"
ref_e = f"{ref_folder}/ref_e.fa"
ref_f = f"{ref_folder}/ref_f.fa"
ref_g = f"{ref_folder}/ref_g.fa"
ref_h = f"{ref_folder}/ref_h.fa"
ref_i = f"{ref_folder}/ref_i.fa"

ref_genomes = {"a": ref_a, "b": ref_b, "c": ref_c, "d": ref_d, "e": ref_e, "f": ref_f, "g": ref_g, "h": ref_h, "i": ref_i}

#print(ref_genomes["a"])

# Location of reads:
reads = ""    #location of reads
#reads = "/Users/xschmd/Desktop/SCOPEWholeGenomeRun/for_analysis/bc17/barcode17_filtered.fastq.gz"