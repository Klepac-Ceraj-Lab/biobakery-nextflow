process metaphlan {
    tag "metaphlan on $sample"
    publishDir "$params.outdir/metaphlan", pattern: "{*.tsv}"

    input:
    tuple val(sample), path(kneads)
    path unmatched
    path metaphlan_db

    output:
    val  sample                  , emit: sample
    path "${sample}_profile.tsv" , emit: profile
    path "${sample}_grouped.fastq.gz"
    path "${sample}_bowtie2.tsv"
    path "${sample}.sam"

    script:
    def forward = kneads[0]
    def reverse = kneads[1]
    def unf = unmatched[0]
    def unr = unmatched[1]

    """
    cat $forward $reverse $unf $unr > ${sample}_grouped.fastq.gz
    
    metaphlan ${sample}_grouped.fastq.gz ${sample}_profile.tsv \
        --bowtie2out ${sample}_bowtie2.tsv \
        --samout ${sample}.sam \
        --input_type fastq \
        --nproc ${task.cpus} \
        --bowtie2db $metaphlan_db
    """
}
 
 process metaphlan_bzip {
    tag "metaphlan_bzip on $sample"
    publishDir "$params.outdir/metaphlan"
    stageInMode "copy"

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