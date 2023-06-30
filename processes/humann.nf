

process humann {
    tag "humann on $sample"
    publishDir "$params.outdir/humann/main"
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }
    cpus { workflow.profile == 'standard' ? null : cpus * task.attempt }

    errorStrategy { task.exitStatus in 134..140 ? 'retry' : 'terminate' }
    maxRetries 3


    input:
    val  sample
    path profile
    path catkneads
    path humann_bowtie_db
    path humann_protein_db

    output:
    val  sample                       , emit: sample
    path "${sample}_genefamilies.tsv" , emit: genefamilies
    path "${sample}_pathabundance.tsv"
    path "${sample}_pathcoverage.tsv"

    script:

    """
    humann_config --update database_folders nucleotide `realpath $humann_bowtie_db`
    humann_config --update database_folders protein `realpath $humann_protein_db`

    humann --input $catkneads --taxonomic-profile $profile --output ./ \
        --threads ${task.cpus} --remove-temp-output --search-mode uniref90 \
        --output-basename $sample
    """
}

process humann_regroup {
    tag "humann_regroup on $sample"
    publishDir "$params.outdir/humann/regroup"

    input:
    val  sample
    path genefamilies
    path humann_utility_db

    output:
    val  sample , emit: sample
    path "${sample}_ecs.tsv"
    path "${sample}_kos.tsv"
    path "${sample}_pfams.tsv"

    script:

    """
    humann_config --update database_folders utility_mapping `realpath $humann_utility_db`
    humann_regroup_table --input $genefamilies --output ${sample}_ecs.tsv --groups uniref90_level4ec
    humann_regroup_table --input $genefamilies --output ${sample}_kos.tsv --groups uniref90_ko
    humann_regroup_table --input $genefamilies --output ${sample}_pfams.tsv --groups uniref90_pfam
    """
}   

process humann_rename {
    tag "humann_rename on $sample"
    publishDir "$params.outdir/humann/rename"

    input:
    val sample
    path ecs
    path kos
    path pfams
    path humann_utility_db

    output:
    val  sample , emit: sample
    path "${sample}_ecs_rename.tsv"
    path "${sample}_kos_rename.tsv"
    path "${sample}_pfams_rename.tsv"

    script:

    """
    humann_config --update database_folders utility_mapping `realpath $humann_utility_db`
    humann_rename_table --input $ecs   --output ${sample}_ecs_rename.tsv   --names ec
    humann_rename_table --input $kos   --output ${sample}_kos_rename.tsv   --names kegg-orthology
    humann_rename_table --input $pfams --output ${sample}_pfams_rename.tsv --names pfam
    """
}