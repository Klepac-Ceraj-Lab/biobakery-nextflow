#!/usr/bin/env nextflow

params.reads = "$baseDir/test/rawfastq/*_L00{1,2,3,4}_R{1,2}_001.fastq.gz"
params.outdir = "$baseDir/output"
params.metaphlan_procs = 8
params.humann_procs = 8

workflow {
    read_pairs_ch = Channel
        .fromFilePairs(params.reads, size: 8)

    knead_out = kneaddata(read_pairs_ch)
    metaphlan_out = metaphlan(knead_out[0], knead_out[1])
    humann(metaphlan_out[0], knead_out[1])
    
}

process kneaddata {
    tag "kneaddata $sample"
    publishDir "$params.outdir/kneaddata"

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample)                                     , emit: sample
    path("${sample}_kneaddata_paired_{1,2}.fastq.gz")     , emit: paired
    path("${sample}_kneaddata*.fastq.gz") , optional:true , emit: others
    path("${sample}_kneaddata.log")                       , emit: log

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
    tuple val(sample)
    path(kneads)

    output:
    tuple val(sample), path("${sample}_profile.tsv") , emit: profile
    path "${sample}_bowtie2.tsv"
    path "${sample}.sam"

    script:
    def forward = kneads[0]
    def reverse = kneads[1]

    """
    cat $forward $reverse > ${sample}_grouped.fastq.gz
    metaphlan ${sample}_grouped.fastq.gz ${sample}_profile.tsv --bowtie2out ${sample}_bowtie2.tsv --samout ${sample}.sam --input_type fastq --nproc ${params.metaphlan_procs}
    """
}
 
process humann {
    publishDir "$params.outdir/humann/main"

    input:
    tuple val(sample), path(profile)
    path kneads

    output:
    tuple val(sample), path("${sample}_genefamilies.tsv") , emit: genefamilies
    path "${sample}_pathabundance.tsv"
    path "${sample}_pathcoverage.tsv"

    script:
    def forward = kneads[0]
    def reverse = kneads[1]

    """
    cat $forward $reverse > ${sample}_grouped.fastq.gz
    humann --input ${sample}_grouped.fastq.gz --taxonomic-profile $profile --output ./ --threads ${params.humann_procs} --remove-temp-output --search-mode uniref90 --output-basename $sample
    """
}