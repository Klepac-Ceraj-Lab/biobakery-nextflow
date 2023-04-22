#!/usr/bin/env nextflow

workflow {
    read_pairs_ch = Channel
        .fromFilePairs("$params.readsdir/$params.filepattern", size: 2)

    knead_out     = kneaddata(read_pairs_ch)
    metaphlan_out = metaphlan(knead_out[0], knead_out[1])
    // humann_out    = humann(metaphlan_out[0], metaphlan_out[1], metaphlan_out[2])
    // regroup_out   = humann_regroup(humann_out[0], humann_out[1])
    // humann_rename(regroup_out)
}

process kneaddata {
    tag "kneaddata $sample"
    publishDir "$params.outdir/kneaddata"
    time { workflow.profile == 'standard' ? null : time * task.attempt }
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }
    
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple val(sample), path(reads)

    output:
    val  sample                                          , emit: sample
    path "${sample}_kneaddata_paired_{1,2}.fastq.gz"     , emit: pairedx
    path "${sample}_kneaddata*.fastq.gz" , optional:true , emit: others
    path "${sample}_kneaddata.log"                       , emit: log

    script:
    
    """
    echo $sample

       kneaddata --input ${reads[0]} --input ${reads[1]} --reference-db /hg37 --output ./ --processes ${task.cpus} --output-prefix ${sample}_kneaddata --trimmomatic /opt/conda/share/trimmomatic

    gzip *.fastq
    """  
}

process metaphlan {
    tag "metaphlan on $sample"
    publishDir "$params.outdir/metaphlan", pattern: "{*.tsv,*.sam}"

    input:
    val sample
    path kneads

    output:
    val  sample                  , emit: sample
    path "${sample}_profile.tsv" , emit: profile
    path "${sample}_grouped.fastq.gz"
    path "${sample}_bowtie2.tsv"
    path "${sample}.sam"

    script:
    def forward = kneads[0]
    def reverse = kneads[1]

    """
    cat $forward $reverse > ${sample}_grouped.fastq.gz
    metaphlan ${sample}_grouped.fastq.gz ${sample}_profile.tsv --bowtie2out ${sample}_bowtie2.tsv --samout ${sample}.sam --input_type fastq --nproc ${task.cpus}
    """
}
 
process humann {
    tag "humann on $sample"
    publishDir "$params.outdir/humann/main"
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }
    cpus { workflow.profile == 'standard' ? null : cpus * task.attempt }
    
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3


    input:
    val  sample
    path profile
    path catkneads

    output:
    val  sample                       , emit: sample
    path "${sample}_genefamilies.tsv" , emit: genefamilies
    path "${sample}_pathabundance.tsv"
    path "${sample}_pathcoverage.tsv"

    script:
    """
    humann --input $catkneads --taxonomic-profile $profile --output ./ --threads ${task.cpus} --remove-temp-output --search-mode uniref90 --output-basename $sample
    """
}

process humann_regroup {
    tag "humann_regroup on $sample"
    publishDir "$params.outdir/humann/regroup"

    input:
    val  sample
    path genefamilies

    output:
    val  sample , emit: sample
    path "${sample}_ecs.tsv"
    path "${sample}_kos.tsv"
    path "${sample}_pfams.tsv"

    script:

    """
    humann_regroup_table --input $genefamilies --output ${sample}_ecs.tsv --groups uniref90_level4ec
    humann_regroup_table --input $genefamilies --output ${sample}_kos.tsv --groups uniref90_ko
    humann_regroup_table --input $genefamilies --output ${sample}_pfams.tsv --groups uniref90_pfam
    """
}   

process humann_rename {
    tag "humann_rename on $sample"
    publishDir "$params.outdir/humann/rename"

    input:
    val  sample
    path ecs
    path kos
    path pfams

    output:
    val  sample , emit: sample
    path "${sample}_ecs_rename.tsv"
    path "${sample}_kos_rename.tsv"
    path "${sample}_pfams_rename.tsv"

    script:

    """
    humann_rename_table --input $ecs   --output ${sample}_ecs_rename.tsv   --names ec
    humann_rename_table --input $kos   --output ${sample}_kos_rename.tsv   --names kegg-orthology
    humann_rename_table --input $pfams --output ${sample}_pfams_rename.tsv --names pfam
    """
}