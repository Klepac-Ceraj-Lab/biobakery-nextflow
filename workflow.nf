#!/usr/bin/env nextflow
 
/*
 * The following pipeline parameters specify the refence genomes
 * and read pairs and can be provided as command line options
 */

params.reads = "$baseDir/test/rawfastq/*_L00{1,2,3,4}_R{1,2}_001.fastq.gz"
params.outdir = "results"


workflow {
    read_pairs_ch = Channel
        .fromFilePairs(params.reads, size: 8)
        .view()
 
    // INDEX(params.transcriptome)
    // FASTQC(read_pairs_ch)
    // QUANT(INDEX.out, read_pairs_ch)
}
 
process INDEX {
    tag "$transcriptome.simpleName"
 
    input:
    path transcriptome
 
    output:
    path 'index'
 
    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}
 
process FASTQC {
    tag "FASTQC on $sample_id"
    publishDir params.outdir
 
    input:
    tuple val(sample_id), path(reads)
 
    output:
    path "fastqc_${sample_id}_logs"
 
    script:
    """
    fastqc.sh "$sample_id" "$reads"
    """
}
 
process QUANT {
    tag "$pair_id"
    publishDir params.outdir
 
    input:
    path index
    tuple val(pair_id), path(reads)
 
    output:
    path pair_id
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}