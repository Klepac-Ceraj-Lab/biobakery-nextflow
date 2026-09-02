process baqlava {
    // Viral profiling using BAQLaVa; runs after humann using kneaddata reads + metaphlan profile
    tag "baqlava on $sample"
    publishDir "$params.outdir/baqlava", mode: 'copy'

    input:
    val(sample)
    path(kneads)
    path(profile)

    output:
    val(sample),                emit: sample
    path("${sample}_baqlava/"), emit: results

    script:
    // --bypass-bacterial-depletion skips the internal humann step.
    // Only needed when the MetaPhlAn prescreen finds 0 species (e.g. tiny test samples);
    // real microbiome samples should keep depletion enabled (default: false).
    def depletion_flag = params.baqlava_bypass_depletion ? "--bypass-bacterial-depletion" : "--taxonomic-profile $profile"

    """
    baqlava -i $kneads \
        $depletion_flag \
        -o ${sample}_baqlava \
        --threads ${task.cpus} \
        --local-jobs 1
    """
}
