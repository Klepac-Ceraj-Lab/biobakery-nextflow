process kneaddata {
    tag "kneaddata $sample"
    publishDir "$params.outdir/kneaddata"
    time { workflow.profile == 'standard' ? null : time * task.attempt }
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }
    maxForks 4

    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(sample), path(reads)
    path human_genome

    output:
    val  sample                                          , emit: sample
    path "${sample}_kneaddata_paired_{1,2}.fastq.gz"     , emit: pairedx
    path "${sample}_kneaddata*.fastq.gz" , optional:true , emit: others
    path "${sample}_kneaddata.log"                       , emit: log

    script:
    
    """
    echo $sample

    kneaddata --input ${reads[0]} --input ${reads[1]} \
              --reference-db $human_genome --output ./ \
              --processes ${task.cpus} --output-prefix ${sample}_kneaddata \
              --trimmomatic /opt/conda/share/trimmomatic

    gzip *.fastq
    """  
}