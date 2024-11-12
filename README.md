# hbv_nano

hbv_nano is a pipeline for analysing whole genome Nanopore long-read data for resistance typing of hepatitis B viruses (HBV).

## Requirements
To run the pipleline `Mamba` and `Apptainer` have to be installed.

## Installation
### 1. Clone the repository and change into it.
Clone the repository
```
git clone -b hbv_nano_daniel https://github.com/ClinicalGenomicsGBG/hbv_nano.git
```
and change the location to it.
```
cd /hbv_nano
```

### 2. Create a mamba base enironment to run the pipeline:
```
mamba env create -f workflow/scripts/environment.yaml
```
If you want to create the environment in a specific location add the flag `-p <path>`
```
mamba env create -f workflow/scripts/environment.yaml -p <path>
```

### 3. Build Apptainer images for the different software needed to run the pipeline.
Make sure that you are still in the location ``.../hbv_nano``. To build the images run the script `build_images.sh`:
```
./build_images.sh
```

### 4. Prepare input fastq.gz files
The pipeline takes gzipped `.fastq.gz` files in the format `{sample_01}.fastq.gz`, `{sample_02}.fastq.gz`, ...,  `{sample_**}.fastq.gz` as input.

The negative control sample should be named like this:
```
{sample}_neg_ctrl.fastq.gz
```

### 5. Set input and ouput
Open the file `config/config.yaml` and set input folder for fastq files under `fastq_folder` and output folder under `output`.
```
# Folder with input .fastq files.
fastq_folder: "data/<folder with input .fastq files>"

# Output folder:
output: "<output folder>"
```

## Running hbv_nano
Depending on your setup and what kind of user you are there are different ways of running the piepline.

### 1. Internal users (those working at BDC)
Run the pipeline using qsub. See internal BDC documentation in Confluence.

### 2. External users
Depending on how your local setup you can run the pipeline in the following way:
#### a) Activate the mamba environment
Activate the environment (matching your naming of the environment):
```
mamba activate hbv_nano
```
##### b) Run the workflow
Change the flag `--cores` according to your resources.

Run the workflow using this command:
```
snakemake --cores 40 --software-deployment-method apptainer --configfile run/config.yaml
```

## Output
Verify that the workflow has finished running sucsessfully by checking that `Finished job 0` appears towards the end of the logfile.

### Locate the output files
When the run is completed all the results can be found in the output folder as set in `config.yaml`. The most important files (depending on your usecase) can be found in the folder `clinic` located in the output folder. The types of files located in the folder `clinic` are:

#### Consensus fasta files
The consenus sequences for the samples in the format `{sample}.ref_*_medaka.fa`.

#### Modified variant calling files
Files in the format `variant_calling_{sample}_ref_*.txt` containing information regarding amino acid changes and which possible drug resistance mutations they may cause. This is listed for all the alternative sequences for each sample.

#### Quality control file
The QC file has the name `qc_clinic.csv`. The file contains the following columns:

| read_id | ref | mapped_reads | mapped_reads_rt | qc_pass |
| ------- | --- | ------------ | --------------- | ------- |

`read_id = the name of the sample, {sample}`

`ref = the reference/genotype which most closely matches the sample`

`mapped_reads = the number of mapped reads for the sample`

`mapped_reads_rt = the number of mapped reads for the sample in the rt region`

`qc_pass = True or False. Whether the sample passed the QC or not.`