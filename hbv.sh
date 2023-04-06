#hbv pipeline (director's cut)_
conda activate /Users/daniel/miniconda/envs/minimap2

ref_a=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/ref_a.fa		#location of reference sequence
ref_b=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/ref_b.fa
ref_c=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/ref_c.fa
ref_d=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/ref_d.fa
ref_e=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/ref_e.fa
ref_f=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/ref_f.fa
ref_g=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/ref_g.fa
ref_h=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/ref_h.fa
ref_i=/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/ref_i.fa

reads=/Users/daniel/Desktop/hbv/sequences_0331_rapid_flongle/KH22-5514_provförväx_kontamination__R1.fastq 	#location of reads

samout_a="${reads:0:-6}_al_a.sam"
samout_b="${reads:0:-6}_al_b.sam"
samout_c="${reads:0:-6}_al_c.sam"
samout_d="${reads:0:-6}_al_d.sam"
samout_e="${reads:0:-6}_al_e.sam"
samout_f="${reads:0:-6}_al_f.sam"
samout_g="${reads:0:-6}_al_g.sam"
samout_h="${reads:0:-6}_al_h.sam"
samout_i="${reads:0:-6}_al_i.sam"

#/Users/daniel/Desktop/hbv/hbv_referensgenom_nya/I_AF241409_Vietnam.

minimap2 -ax map-ont $ref_a $reads > $samout_a
minimap2 -ax map-ont $ref_b $reads > $samout_b
minimap2 -ax map-ont $ref_c $reads > $samout_c
minimap2 -ax map-ont $ref_d $reads > $samout_d
minimap2 -ax map-ont $ref_e $reads > $samout_e
minimap2 -ax map-ont $ref_f $reads > $samout_f
minimap2 -ax map-ont $ref_g $reads > $samout_g
minimap2 -ax map-ont $ref_h $reads > $samout_h
minimap2 -ax map-ont $ref_i $reads > $samout_i

conda deactivate

conda activate /Users/daniel/miniconda/envs/samtools

sort_a="${samout_a:0:-4}_sorted.bam"
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

 samtools mpileup
