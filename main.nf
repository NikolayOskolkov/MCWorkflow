#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define absolute paths to pseudo-reads and annotation
def pseudo_reads_file = file(params.input_pseudo_reads)
def output_dir        = file(params.input_dir)

Channel
    .fromPath("${params.input_dir}/*")
    .map { file ->
        tuple(
            file,
            output_dir,
            params.type_of_pseudo_reads,
            params.threads,
            pseudo_reads_file,
            params.n_allowed_multimappers
        )
    }
    .set { input_data }

workflow {
    indexed = index_reference(input_data)
    aligned = align_pseudo_reads(indexed)
    detect_exogenous(aligned)
}

// Process 1: Indexing
process index_reference {
    input:
    tuple path(input_ref), path(output_dir), val(type_of_pseudo_reads), val(threads), path(input_pseudo_reads), val(n_allowed_multimappers)

    output:
    tuple path(input_ref), path(output_dir), val(type_of_pseudo_reads), val(threads), path(input_pseudo_reads), val(n_allowed_multimappers)

    script:
    """
    cd "${output_dir}"
    bowtie2-build --large-index ${input_ref} ${input_ref} --threads ${threads} >> bowtie2-build.log 2>&1
    """
}

// Process 2: Alignment
process align_pseudo_reads {
    input:
    tuple path(input_ref), path(output_dir), val(type_of_pseudo_reads), val(threads), path(input_pseudo_reads), val(n_allowed_multimappers)

    output:
    tuple path(input_ref), path(output_dir), val(type_of_pseudo_reads), val(threads), path(input_pseudo_reads), val(n_allowed_multimappers)

    script:
    """
    cd "${output_dir}"
    mkdir -p ${input_ref}_${type_of_pseudo_reads}
    cd ${input_ref}_${type_of_pseudo_reads}

    bowtie2 --large-index -f -k ${n_allowed_multimappers} -x ../${input_ref} \
        --end-to-end --quiet --threads ${threads} --very-sensitive \
        -U ../../${input_pseudo_reads} | \
        samtools view -bS -F 4 -h -@ ${threads} - | \
        samtools sort -@ ${threads} - > PseudoReads_aligned_to_${input_ref}.bam

    samtools index -c PseudoReads_aligned_to_${input_ref}.bam
    """
}

// Process 3: Detection
process detect_exogenous {
  input:
    tuple path(input_ref), path(output_dir), val(type_of_pseudo_reads), val(threads), path(input_pseudo_reads), val(n_allowed_multimappers)

  script:
  """
  bash ${projectDir}/detect_exogenous.sh \
      ${input_ref} \
      ${output_dir} \
      ${type_of_pseudo_reads}
  """
}

