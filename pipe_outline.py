"""

#1. fastp to trim reads, remove adapters and filter out low quality reads
fastp -i <reads.fastq> -o <reads_trimmed.fastq> -h <fastp.html>

#2. (done in hbv.sh)
minimap2 -ax map-ont ref.fasta reads.fastq > alignment.sam  # map reads to references (a,b,c,d,e,f,g,h,i)
                                                            # this produces 9 different alignments
#3. (done in hbv.sh)
samtools view -S -b alignment.sam | samtools sort -o sorted_alignment.bam   # convert .sam files to .bam and sort them

#4. (done in hbv.sh)
samtools stats sorted_alignment.bam |grep -i 'error rate' | cut -f 3    # get error rate (used for selecting best alignment). Take the lowest error rate.


#5. - generate consensus sequence (cluster)
samtools consensus -o <consensus.fa> -m simple -d20 <sorted_alignment.bam>    # generate consensus sequence from sorted alignment

#6. - Polish consensus sequence (cluster)
medaka_consensus -i <reads.fastq.gz> -d <consensus.fa> -o <output folder>  -m r941_min_sup_g507    # polish consensus
r941_min_sup_g507 = {pore}_{device}_{caller variant}_{caller version}

    qlogin -q development.q -pe mpi 5

    module load miniconda/4.14.0

    source activate /home/xschmd/.conda/envs/medaka_py_3.8

    medaka_consensus -i <reads.fastq.gz> -d <consensus.fa> -o <output folder> -m r1041_e82_400bps_sup_g615

"""