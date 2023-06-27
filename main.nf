#!/usr/bin/env nextflow
params.assmblreads = '/Users/daniel/Desktop/hbv/validering/hbv_val_05/barcode01/bc01.fastq.gz'
params.assmblref = '/Users/daniel/Desktop/hbv/hbv_referensgenom/ref_a.fa'


process minimap2 {
    conda 'bioconda::minimap2'    
    input:
    tuple val(readID), path(readFile)
    path genome

    output:
    tuple val(readID), path("${readID}.paf")

    """
    minimap2 "${readFile}" "${genome}" > "${readID}.paf"
    """
}

workflow {

    refgenome_file = file( params.assmblref )
    fasta_files = Channel.fromPath( params.assmblreads )
                         .map { file -> tuple( file.baseName, file ) }

    minimap2( fasta_files, refgenome_file )
}