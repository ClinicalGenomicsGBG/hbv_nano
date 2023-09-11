# This is deprecated and will be replaced by samtools consensus

# Use source consensus.sh to run script
# Run this script after hbv.sh has been run

conda activate /Users/xschmd/miniconda3/hbv

ref_a=/Users/xschmd/Desktop/referensgenom/ref_a.fa  #location of reference sequences
ref_b=/Users/xschmd/Desktop/referensgenom/ref_b.fa
ref_c=/Users/xschmd/Desktop/referensgenom/ref_c.fa
ref_d=/Users/xschmd/Desktop/referensgenom/ref_d.fa
ref_e=/Users/xschmd/Desktop/referensgenom/ref_e.fa
ref_f=/Users/xschmd/Desktop/referensgenom/ref_f.fa
ref_g=/Users/xschmd/Desktop/referensgenom/ref_g.fa
ref_h=/Users/xschmd/Desktop/referensgenom/ref_h.fa
ref_i=/Users/xschmd/Desktop/referensgenom/ref_i.fa

#sorted=/Users/daniel/Desktop/hbv/validering/hbv_val_05/barcode06/bc06.fa_al_d_sorted.bam     #sorted bam file with reads
sorted=/Users/xschmd/Desktop/validering/hbv_val_07/fastq_pass/barcode13/bc13_filtered.fa_al_c_sorted.bam
out="${sorted:0:-11}_consensus_d20.fasta"

samtools mpileup -uf  $ref_c $sorted | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq -d20 > $out     #change ref_x   #can add -d20 for vcfutils to only include reads with depth >20
                                                                                                        #standard is -d3
#bcftools mpileup -f $ref_d $sorted | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq -d20 > $out     #change ref_x   #can add -d20 for vcfutils to only include reads with depth >20


#bcftools mpileup -Ou -f $ref_c $sorted | bcftools call -mv -Oz -o /Users/xschmd/Desktop/test2.vcf.gz
#bcftools index /Users/xschmd/Desktop/test2.vcf.gz
#bcftools consensus -f $ref_c /Users/xschmd/Desktop/test2.vcf.gz > /Users/xschmd/Desktop/test2.fa

echo $out

conda deactivate


