#must use source consensus.sh to run script

#run this script after hbv.sh has been run
conda activate /Users/daniel/miniconda/envs/samtools

ref_a=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_a.fa		#location of reference sequence
ref_b=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_b.fa
ref_c=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_c.fa
ref_d=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_d.fa
ref_e=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_e.fa
ref_f=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_f.fa
ref_g=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_g.fa
ref_h=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_h.fa
ref_i=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_i.fa

sorted=/Users/daniel/Desktop/hbv/validering/hbv_val_04/barcode23/bc23_fastp.fa_al_e_sorted.bam        #sorted bam file with reads
out="${sorted:0:-11}_consensus.fasta"

samtools mpileup -uf  $ref_e $sorted | bcftools call -c --ploidy 1 | vcfutils.pl vcf2fq > $out     #change ref_x

echo $out

conda deactivate


