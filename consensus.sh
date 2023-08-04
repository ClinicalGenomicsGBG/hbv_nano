#must use source consensus.sh to run script

#run this script after hbv.sh has been run
#conda activate /Users/daniel/miniconda/envs/samtools

conda activate /Users/xschmd/miniconda3/hbv
#ref_a=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_a.fa		#location of reference sequence
#ref_b=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_b.fa
#ref_c=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_c.fa
#ref_d=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_d.fa
#ref_e=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_e.fa
#ref_f=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_f.fa
#ref_g=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_g.fa
#ref_h=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_h.fa
#ref_i=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_i.fa

ref_a=/Users/xschmd/Desktop/referensgenom/ref_a.fa
ref_b=/Users/xschmd/Desktop/referensgenom/ref_b.fa
ref_c=/Users/xschmd/Desktop/referensgenom/ref_c.fa
ref_d=/Users/xschmd/Desktop/referensgenom/ref_d.fa
ref_e=/Users/xschmd/Desktop/referensgenom/ref_e.fa
ref_f=/Users/xschmd/Desktop/referensgenom/ref_f.fa
ref_g=/Users/xschmd/Desktop/referensgenom/ref_g.fa
ref_h=/Users/xschmd/Desktop/referensgenom/ref_h.fa
ref_i=/Users/xschmd/Desktop/referensgenom/ref_i.fa

#sorted=/Users/daniel/Desktop/hbv/validering/hbv_val_05/barcode06/bc06.fa_al_d_sorted.bam     #sorted bam file with reads
sorted=/Users/xschmd/Desktop/validering/hbv_val_06/fastq_pass/barcode11/bc11_filtered.fa_al_d_sorted.bam
out="${sorted:0:-11}_consensus_d20.fasta"

#samtools mpileup -uf  $ref_d $sorted | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq -d1000000000 > $out     #change ref_x   #can add -d20 for vcfutils to only include reads with depth >20
                                                                                                        #standard is -d3
bcftools mpileup -f $ref_d $sorted | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq -d20 > $out     #change ref_x   #can add -d20 for vcfutils to only include reads with depth >20


#bcftools mpileup -Ou -f $ref_d $sorted | bcftools call -mv -Oz -o /Users/xschmd/Desktop/test.vcf.gz
#bcftools index /Users/xschmd/Desktop/test.vcf.gz
#bcftools consensus -f $ref_d /Users/xschmd/Desktop/test.vcf.gz > /Users/xschmd/Desktop/test.fa

echo $out

conda deactivate


