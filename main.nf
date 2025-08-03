#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_MCWorkflow {

    input:
    path input_ref
    path input_dir
    val type_of_pseudo_reads
    val threads
    path input_pseudo_reads
    path input_gtdb_annot

    script:
    """
    bash ${projectDir}/micr_cont_detect.sh $input_ref $input_dir $type_of_pseudo_reads $threads $input_pseudo_reads $input_gtdb_annot
    """
}

workflow {
    run_MCWorkflow(file(params.input_ref), file(params.input_dir), params.type_of_pseudo_reads, params.threads, file(params.input_pseudo_reads), file(params.input_gtdb_annot))
}

