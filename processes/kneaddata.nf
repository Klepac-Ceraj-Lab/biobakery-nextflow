process single_end_kneaddata {
    // Run kneaddata on single reads
    tag "kneaddata $sample"
    publishDir "$params.outdir/kneaddata", mode: 'link'
    time { workflow.profile == 'standard' ? null : time * task.attempt }
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }

    when:
    params.paired_end == false

    input:
    tuple val(sample), path(reads)

    output:
    val(sample), emit: sample
    path("${sample}_kneaddata.fastq.gz"), emit: kneads
    path "${sample}_kneaddata*.fastq.gz" , optional:true , emit: others
    path "${sample}_kneaddata.log"                       , emit: log

    shell:
    
    """
    echo $sample

    kneaddata --unpaired $reads \
              --reference-db ${params.host_genome} --output ./ \
              --threads ${task.cpus} --output-prefix ${sample}_kneaddata

    gzip ${sample}_kneaddata*.fastq
    """  
}

process paired_end_kneaddata {
    // Run kneaddata on paired reads
    tag "kneaddata $sample"
    publishDir "$params.outdir/kneaddata", mode: 'link', saveAs: { f ->
      // filenames to exclude from publishing
      def n = f.toString()
      if ( n ==~ /.+_concatenated\.fastq\.gz/ ) return null
      // publish everything else
      return n
    }
    time { workflow.profile == 'standard' ? null : time * task.attempt }
    memory { workflow.profile == 'standard' ? null : memory * task.attempt }

    when:
    params.paired_end == true

    input:
    tuple val(sample), path(reads)

    output:
    val(sample), emit: sample
    path("${sample}_kneaddata_paired_1.fastq.gz"), emit: paired1
    path("${sample}_kneaddata_paired_2.fastq.gz"), emit: paired2
    path("${sample}_kneaddata_unmatched_1.fastq.gz"), emit: unpaired1
    path("${sample}_kneaddata_unmatched_2.fastq.gz"), emit: unpaired2
    path "${sample}_kneaddata.log"                       , emit: log
    path "${sample}_concatenated.fastq.gz"                       , emit: kneads // available for chaining; not published

    shell:
    
    """

    echo $sample
    echo $reads

    kneaddata -i1 ${reads[0]} -i2 ${reads[1]} \
              --reference-db ${params.host_genome} --output ./ \
              --processes ${task.cpus} --output-prefix ${sample}_kneaddata

    gzip ${sample}_kneaddata*.fastq
    
    # concatenate singleton/paired reads into a file to feed into metaphlan
    cat ${sample}_kneaddata_paired*.fastq.gz ${sample}_kneaddata_unmatched*.fastq.gz > ${sample}_concatenated.fastq.gz
    """



}
