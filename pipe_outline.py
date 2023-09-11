"""

#0.5. fastp to trim reads, remove adapters and filter out low quality reads
fastp -i reads.fastq -o reads_trimmed.fastq -h fastp.html

#1. (done in hbv.sh)
minimap2 -ax map-ont ref.fasta reads.fastq > alignment.sam  # map reads to references (a,b,c,d,e,f,g,h,i)
                                                            # this produces 9 different alignments
#2.
samtools view -S -b alignment.sam | samtools sort -o sorted_alignment.bam   # convert .sam files to .bam and sort them

#3.
samtools stats sorted_alignment.bam |grep -i 'error rate' | cut -f 3    # get error rate (used for selecting best alignment). Take the lowest error rate.

#4. - generate consensus sequence
samtools consensus -o <consensus.fa> -m simple -d20 <sorted_alignment.bam>    # generate consensus sequence from sorted alignment

#5.
medaka_consensus -i reads.fastq -d consensus.fasta -o /outpath  -m r941_min_sup_g507    # use medaka to polish consensus
r941_min_sup_g507 = {pore}_{device}_{caller variant}_{caller version}

"""
