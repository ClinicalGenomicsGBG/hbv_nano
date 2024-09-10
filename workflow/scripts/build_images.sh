#!/bin/bash -l

apptainer build images/fastp_0.23.4.img docker://staphb/fastp
apptainer build images/freebayes_1.3.7.img docker://staphb/freebayes
apptainer build images/medaka_1.12.0.img docker://ontresearch/medaka
apptainer build images/minimap2_2.28.img docker://staphb/minimap2
apptainer build images/samtools_1.20.img docker://staphb/samtools