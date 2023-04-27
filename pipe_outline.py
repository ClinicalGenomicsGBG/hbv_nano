"""

#1.
minimap2 -ax map-ont ref.fasta reads.fastq > alignment.sam  # map reads to reference

#2.
samtools view -S -b alignment.sam | samtools sort -o sorted_alignment.bam   # convert to bam and sort

#3.
samtools stats sorted_alignment.bam |grep -i 'error rate' | cut -f 3    # get error rate (used for selecting best alignment)

#4.
samtools mpileup -uf ref_genome.fasta sorted_alignment.bam | bcftools call -c | vcfutils.pl vcf2fq > consensus.fasta # get consensus

#5.
medaka_consensus -i reads.fastq -d consensus.fasta -o /outpath  -m r941_min_sup_g507 # use medaka to polish consensus
r941_min_sup_g507 = {pore}_{device}_{caller variant}_{caller version}

"""