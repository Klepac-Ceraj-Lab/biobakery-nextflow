#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { metaphlan; metaphlan_bzip } from './processes/metaphlan.nf'
include { humann; humann_regroup; humann_rename } from './processes/humann.nf'

workflow {
    
    read_pairs_ch = Channel
        .fromFilePairs("$params.readsdir/$params.filepattern", size: 2)

    human_genome      = params.human_genome
    metaphlan_db      = params.metaphlan_db 
    humann_bowtie_db  = params.humann_bowtie_db 
    humann_protein_db = params.humann_protein_db 
    humann_utility_db = params.humann_utility_db 
    
    knead_out     = kneaddata(read_pairs_ch, human_genome)
    metaphlan_out = metaphlan(knead_out[0], metaphlan_db)
    metaphlan_bzip = metaphlan_bzip(metaphlan_out[0], metaphlan_out[4])
    humann_out    = humann(metaphlan_out[0], metaphlan_out[1], metaphlan_out[2], humann_bowtie_db, humann_protein_db)
    regroup_out   = humann_regroup(humann_out[0], humann_out[1], humann_utility_db)
    humann_rename(regroup_out, humann_utility_db)
}

process kneaddata {
    tag "kneaddata $sample"
    publishDir "$params.outdir/kneaddata"
    time { workflow.profile == 'standard' ? null : time * task.attempt }
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }

    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample), path(reads)
    path human_genome

    output:
    tuple val(sample) path("${sample}_kneaddata_paired_{1,2}.fastq.gz")
    path "${sample}_kneaddata*.fastq.gz" , optional:true , emit: others
    path "${sample}_kneaddata.log"                       , emit: log

    script:
    
    """
    echo $sample

    kneaddata --input1 ${reads[0]} --input2 ${reads[1]} \
              --reference-db $human_genome --output ./ \
              --processes ${task.cpus} --output-prefix ${sample}_kneaddata \
              --trimmomatic /opt/conda/bin

    gzip *.fastq
    """  
}