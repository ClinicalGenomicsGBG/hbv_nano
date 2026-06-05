#!/bin/bash -l

apptainer build images/fastp_1.3.3.img docker://staphb/fastp:1.3.3
apptainer build images/freebayes_1.3.7.img docker://staphb/freebayes:1.3.7
apptainer build images/medaka_1.12.0.img docker://ontresearch/medaka:1.12.0
apptainer build images/minimap2_2.28.img docker://staphb/minimap2:2.28
apptainer build images/samtools_1.20.img docker://staphb/samtools:1.20
