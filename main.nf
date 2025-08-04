#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define absolute paths to pseudo-reads and annotation
def pseudo_reads_file = file(params.input_pseudo_reads)
def output_dir        = file(params.input_dir)  // this is the "data/" folder

// Create a channel of all files in the input dir
Channel
    .fromPath("${params.input_dir}/*")
    .map { file ->
        tuple(
            file,                          // input_ref
            output_dir,                    // absolute path to "data/"
            params.type_of_pseudo_reads,
            params.threads,
            pseudo_reads_file,
            params.n_allowed_multimappers
        )
    }
    .set { input_data }

workflow {
    run_MCWorkflow(input_data)
}

process run_MCWorkflow {

    input:
    tuple path(input_ref), path(output_dir), val(type_of_pseudo_reads), val(threads), path(input_pseudo_reads), val(n_allowed_multimappers)

    script:
    """
    # Write output files directly into the "data" folder
    bash ${projectDir}/micr_cont_detect.sh \
        $input_ref \
        $output_dir \
        $type_of_pseudo_reads \
        $threads \
        $input_pseudo_reads \
        $n_allowed_multimappers
    """
}

