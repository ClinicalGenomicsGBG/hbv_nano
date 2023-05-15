#must use source consensus.sh to run script

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

sorted=/Users/daniel/Desktop/hbv/validering/hbv_val_02/val_02/20230511_1357_MN29974_AOK137_0a7c36c2/fastq_pass/barcode12/neg_ctrl_bc12.fa_al_d_sorted.bam
out="${sorted:0:-11}_consensus.fasta"

samtools mpileup -uf  $ref_d $sorted | bcftools call -c | vcfutils.pl vcf2fq > $out

echo $out

conda deactivate


