#!/usr/bin/env nextflow

params.reads = "$baseDir/test/rawfastq/*_L00{1,2,3,4}_R{1,2}_001.fastq.gz"
params.outdir = "$baseDir/output"
params.procs = 8

workflow {
    read_pairs_ch = Channel
        .fromFilePairs(params.reads, size: 8)

    knead_out = kneaddata(read_pairs_ch)
    metaphlan(knead_out[0])
    
    // metaphlan(knead_pairs)
}

process kneaddata {
    tag "kneaddata $sample"
    publishDir "$params.outdir/kneaddata"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_kneaddata_paired_{1,2}.fastq.gz") , emit: paired
    path("${sample}_kneaddata*.fastq.gz")                , optional:true , emit: others
    path("${sample}_kneaddata.log")                                      , emit: log

    script:
    def forward = reads.findAll{ read-> read =~ /.+_R1_.+/ }
    def reverse = reads.findAll{ read-> read =~ /.+_R2_.+/ }
    
    """
    echo $sample

    cat ${forward.join(" ")} > ${sample}_1.fastq.gz
    cat ${reverse.join(" ")} > ${sample}_2.fastq.gz

    kneaddata --input ${sample}_1.fastq.gz --input ${sample}_2.fastq.gz --reference-db /hg37 --output ./ --output-prefix ${sample}_kneaddata --trimmomatic /opt/conda/share/trimmomatic

    gzip *.fastq
    """  
}

process metaphlan {
    // tag = "metaphlan on $sample"
    publishDir "$params.outdir/metaphlan"

    input:
    tuple val(sample), path(kneads)

    output:
    tuple val(sample), path("${sample}_profile.tsv") , emit: profile
    path "${sample}_bowtie2.tsv"
    path "${sample}.sam"

    script:
    def forward = kneads[0]
    def reverse = kneads[1]

    """
    cat $forward $reverse > ${sample}_grouped.fastq.gz
    metaphlan ${sample}_grouped.fastq.gz ${sample}_profile.tsv --bowtie2out ${sample}_bowtie2.tsv --samout ${sample}.sam --input_type fastq --nproc ${params.procs}
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