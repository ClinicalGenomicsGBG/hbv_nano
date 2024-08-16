#!/bin/bash -l

# Exit on any error
set -e

# Define path to images
script_path=`dirname $(realpath $0)`
root_image_path="${script_path}/../workflow/images"
mkdir -p "${root_image_path}"

# Run image build commands
apptainer build "${root_image_path}/fastp_0.23.4.img" docker://staphb/fastp:0.23.4
apptainer build "${root_image_path}/freebayes_1.3.7.img" docker://staphb/freebayes:1.3.7
apptainer build "${root_image_path}/medaka_1.12.0.img" docker://ontresearch/medaka:1.12.0
apptainer build "${root_image_path}/minimap2_2.28.img" docker://staphb/minimap2:2.28
apptainer build "${root_image_path}/samtools_1.20.img" docker://staphb/samtools:1.20
