process INDEX_REFERENCE {

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    tuple val(ID), path(input_ref)

    output:
    tuple val(ID), path("*.bt2l")

    script:
    """
    bowtie2-build --large-index \$(basename ${input_ref}) \$(basename ${input_ref}) --threads "${task.cpus}"
    """
}
