#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { metaphlan; metaphlan_bzip } from './processes/metaphlan.nf'
include { humann; humann_regroup; humann_rename } from './processes/humann.nf'

params.filepattern = "*_kneaddata_paired_{1,2}.fastq.gz"

workflow {
    
    read_pairs_ch = Channel
        .fromFilePairs("$params.readsdir/$params.filepattern", size: 2)

    metaphlan_db      =  params.metaphlan_db 
    humann_bowtie_db  =  params.humann_bowtie_db 
    humann_protein_db =  params.humann_protein_db 
    humann_utility_db =  params.humann_utility_db 
    
    metaphlan_out = metaphlan(read_pairs_ch, metaphlan_db)
    // metaphlan_bzip = metaphlan_bzip(metaphlan_out[0], metaphlan_out[4])
    humann_out    = humann(metaphlan_out[0], metaphlan_out[1], metaphlan_out[2], humann_bowtie_db, humann_protein_db)
    regroup_out   = humann_regroup(humann_out[0], humann_out[1], humann_utility_db)
    humann_rename(regroup_out, humann_utility_db)
}
