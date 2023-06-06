#!/usr/bin/env nextflow

nextflow.enable.dsl=2

workflow {
    knead_out = Channel
        .fromFilePairs("$params.readsdir/$params.filepattern", size: 2)

    metaphlan_out = metaphlan(knead_out)
    metaphlan_bzip = metaphlan_bzip(metaphlan_out[0], metaphlan_out[4])
    humann_out    = humann(metaphlan_out[0], metaphlan_out[1], metaphlan_out[2])
    regroup_out   = humann_regroup(humann_out[0], humann_out[1])
    humann_rename(regroup_out)
}


process metaphlan {
    tag "metaphlan on $sample"
    publishDir "$params.outdir/metaphlan", pattern: "{*.tsv,*.sam}"
    maxForks 2

    input:
    tuple val(sample), path(kneads)

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
 
 process metaphlan_bzip {
    tag "metaphlan_bzip on $sample"
    publishDir "$params.outdir/metaphlan"
    maxForks 2

    input:
    val sample
    path sam

    output:
    val  sample                  , emit: sample
    path "${sample}.sam.bz2"

    script:
    """
    bzip2 -v $sam
    """
}

process humann {
    tag "humann on $sample"
    publishDir "$params.outdir/humann/main"
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }
    cpus { workflow.profile == 'standard' ? null : cpus * task.attempt }
    maxForks 4

    errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
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
    // --nucleotide-database $choco --protein-database $protein
}

process humann_regroup {
    tag "humann_regroup on $sample"
    publishDir "$params.outdir/humann/regroup"
    maxForks 2

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
    maxForks 2

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