#hbv pipeline (director's cut)
#use source hbv.sh to run script
#run consensus.sh after this script

conda activate /Users/xschmd/miniconda3/hbv

ref_a=/Users/xschmd/Desktop/referensgenom/ref_a.fa    #location of reference sequence
ref_b=/Users/xschmd/Desktop/referensgenom/ref_b.fa
ref_c=/Users/xschmd/Desktop/referensgenom/ref_c.fa
ref_d=/Users/xschmd/Desktop/referensgenom/ref_d.fa
ref_e=/Users/xschmd/Desktop/referensgenom/ref_e.fa
ref_f=/Users/xschmd/Desktop/referensgenom/ref_f.fa
ref_g=/Users/xschmd/Desktop/referensgenom/ref_g.fa
ref_h=/Users/xschmd/Desktop/referensgenom/ref_h.fa
ref_i=/Users/xschmd/Desktop/referensgenom/ref_i.fa

reads=/Users/xschmd/Desktop/validering/hbv_val_07/fastq_pass/barcode18/bc18_filtered.fastq.gz

samout_a="${reads:0:-6}_al_a.sam"   #writes the .sam output to the same path as the read input
samout_b="${reads:0:-6}_al_b.sam"
samout_c="${reads:0:-6}_al_c.sam"
samout_d="${reads:0:-6}_al_d.sam"
samout_e="${reads:0:-6}_al_e.sam"
samout_f="${reads:0:-6}_al_f.sam"
samout_g="${reads:0:-6}_al_g.sam"
samout_h="${reads:0:-6}_al_h.sam"
samout_i="${reads:0:-6}_al_i.sam"


minimap2 -ax map-ont $ref_a $reads > $samout_a    #map the reads to the reference sequences (of the different genotypes)
minimap2 -ax map-ont $ref_b $reads > $samout_b
minimap2 -ax map-ont $ref_c $reads > $samout_c
minimap2 -ax map-ont $ref_d $reads > $samout_d
minimap2 -ax map-ont $ref_e $reads > $samout_e
minimap2 -ax map-ont $ref_f $reads > $samout_f
minimap2 -ax map-ont $ref_g $reads > $samout_g
minimap2 -ax map-ont $ref_h $reads > $samout_h
minimap2 -ax map-ont $ref_i $reads > $samout_i

#conda deactivate

#conda activate /Users/daniel/miniconda/envs/samtools
sort_a="${samout_a:0:-4}_sorted.bam"   #naming of the sorted bam files
sort_b="${samout_b:0:-4}_sorted.bam"
sort_c="${samout_c:0:-4}_sorted.bam"
sort_d="${samout_d:0:-4}_sorted.bam"
sort_e="${samout_e:0:-4}_sorted.bam"
sort_f="${samout_f:0:-4}_sorted.bam"
sort_g="${samout_g:0:-4}_sorted.bam"
sort_h="${samout_h:0:-4}_sorted.bam"
sort_i="${samout_i:0:-4}_sorted.bam"

samtools view -S -b $samout_a | samtools sort -o $sort_a
samtools view -S -b $samout_b | samtools sort -o $sort_b
samtools view -S -b $samout_c | samtools sort -o $sort_c
samtools view -S -b $samout_d | samtools sort -o $sort_d
samtools view -S -b $samout_e | samtools sort -o $sort_e
samtools view -S -b $samout_f | samtools sort -o $sort_f
samtools view -S -b $samout_g | samtools sort -o $sort_g
samtools view -S -b $samout_h | samtools sort -o $sort_h
samtools view -S -b $samout_i | samtools sort -o $sort_i

echo "a:"
samtools stats $sort_a |grep -i 'error rate' | cut -f 3
echo "b:"
samtools stats $sort_b |grep -i 'error rate' | cut -f 3
echo "c:"
samtools stats $sort_c |grep -i 'error rate' | cut -f 3
echo "d:"
samtools stats $sort_d |grep -i 'error rate' | cut -f 3
echo "e:"
samtools stats $sort_e |grep -i 'error rate' | cut -f 3
echo "f:"
samtools stats $sort_f |grep -i 'error rate' | cut -f 3
echo "g:"
samtools stats $sort_g |grep -i 'error rate' | cut -f 3
echo "h:"
samtools stats $sort_h |grep -i 'error rate' | cut -f 3
echo "i:"
samtools stats $sort_i |grep -i 'error rate' | cut -f 3

#samtools mpileup
#samtools mpileup -uf /Users/daniel/Desktop/hbv/hbv_referensgenom/ref_e.fa /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode01/bc01.fa_al_e_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode01/bc01.fa_al_e_sorted_consensus.fasta
#samtools mpileup -uf /Users/daniel/Desktop/hbv/hbv_referensgenom/ref_a.fa /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode02/bc02.fa_al_a_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode02/bc02.fa_al_a_sorted_consensus.fasta
#samtools mpileup -uf /Users/daniel/Desktop/hbv/hbv_referensgenom/ref_a.fa /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode03/bc03.fa_al_a_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode03/bc03.fa_al_a_sorted_consensus.fasta
#samtools mpileup -uf /Users/daniel/Desktop/hbv/hbv_referensgenom/ref_d.fa /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode04/bc04.fa_al_d_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode04/bc04.fa_al_d_sorted_consensus.fasta
#samtools mpileup -uf /Users/daniel/Desktop/hbv/hbv_referensgenom/ref_d.fa /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode05/bc05.fa_al_d_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode05/bc05.fa_al_d_sorted_consensus.fasta
#samtools mpileup -uf /Users/daniel/Desktop/hbv/hbv_referensgenom/ref_e.fa /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode06/bc06.fa_al_e_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > /Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode06/bc06.fa_al_e_sorted_consensus.fasta   