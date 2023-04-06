conda activate /Users/daniel/miniconda/envs/samtools

ref_a=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/A2_L13994.fa		#location of reference sequence
ref_b=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/B2_AF121251_Vietnam.fa
ref_c=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/C2_X04615_Japan.fa
ref_d=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/D1_AF121241_Turkey.fa
ref_e=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/E_AB091256.fa
ref_f=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/F_X75663.fa
ref_g=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/G__AB064312.fa
ref_h=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/H_AY090454.fa
ref_i=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/I_AF241409_Vietnam.fa

sorted=/Users/daniel/Desktop/hbv/sequences_0331_rapid_flongle/selected/KH22-5514_provförväx_kontamination__R1_al_a_sorted.bam
out="${sorted:0:-11}_consensus.fasta"

samtools mpileup -uf  $ref_a $sorted | bcftools call -c | vcfutils.pl vcf2fq > $out

echo $out

conda deactivate


