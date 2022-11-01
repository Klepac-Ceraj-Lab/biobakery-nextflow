#!/usr/bin/env nextflow

params.reads = "$baseDir/test/rawfastq/*_L00{1,2,3,4}_R{1,2}_001.fastq.gz"
params.outdir = "$baseDir/output"

workflow {
    read_pairs_ch = Channel
        .fromFilePairs(params.reads, size: 8)

    kneaddata(read_pairs_ch) | view
}

process kneaddata {
    tag "$sample"
    publishDir params.outdir

    input:
    tuple val(sample), path(reads)

    output:
    stdout

    script:
    def forward = reads.findAll{ read-> read =~ /.+_R1_.+/ }
    def reverse = reads.findAll{ read-> read =~ /.+_R2_.+/ }
    
    """
    echo $sample

    cat ${forward.join(" ")} > ${sample}_1.fastq.gz
    cat ${reverse.join(" ")} > ${sample}_2.fastq.gz

    """  
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