#!/bin/bash

#$ -N consensus
#$ -q development.q
#$ -cwd
#$ -o test.out -e test.err
#$ -M daniel.schmidt@gu.se
#$ -m bea
#$ -pe mpi 40


module load samtools
samtools consensus -m simple -d20 <reads.bam> -o <consensus.fa>
