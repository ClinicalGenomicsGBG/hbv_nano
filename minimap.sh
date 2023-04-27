conda activate /Users/daniel/miniconda/envs/minimap2

ref_a=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_a.fa		#location of reference sequence
ref_b=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_b.fa
ref_c=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_c.fa
ref_d=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_d.fa
ref_e=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_e.fa
ref_f=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_f.fa
ref_g=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_g.fa
ref_h=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_h.fa
ref_i=/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_i.fa

reads=/Users/daniel/Desktop/hbv_val_01/1-2-3-4-5-ctrl/20230420_1258_MN29974_AOJ936_7461b1af/fastq_pass/barcode01/bc01.fastq.gz

samout_a="${reads:0:-6}_al_a.sam"   #writes the .sam output to the same path as the read input
samout_b="${reads:0:-6}_al_b.sam"
samout_c="${reads:0:-6}_al_c.sam"
samout_d="${reads:0:-6}_al_d.sam"
samout_e="${reads:0:-6}_al_e.sam"
samout_f="${reads:0:-6}_al_f.sam"
samout_g="${reads:0:-6}_al_g.sam"
samout_h="${reads:0:-6}_al_h.sam"
samout_i="${reads:0:-6}_al_i.sam"

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