process ALIGN_PSEUDO_READS {

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    tuple val(ID), path(index), path(input_pseudo_reads), val(type_of_pseudo_reads), val(n_allowed_multimappers)

    output:
    tuple val(ID), path("*.bam")

    script:
    """
    index1=\$(printf '%s\n' *.bt2* | head -n1)
    ref_name=\$(echo \$index1 | sed 's/.1.bt2l//' | sed 's/.1.bt2//')

    input_pseudo_reads_name=\$(basename "$input_pseudo_reads" | sed -E 's/\\.(fna|fa|fasta)(\\.gz)?\$//')

    bowtie2 --large-index -f -k ${n_allowed_multimappers} -x \${ref_name} \
        --end-to-end --quiet --threads "${task.cpus}" --very-sensitive \
        -U ${input_pseudo_reads} | \
        samtools view -bS -F 4 -h -@ "${task.cpus}" - | \
        samtools sort -@ "${task.cpus}" - > PseudoReads_aligned_to_\${input_pseudo_reads_name}.bam
    
    input_ref_name=\$(echo \$index1 | sed 's/.1.bt2l//' | sed -E 's/.fasta|.fa|.fna//' | sed 's/.gz//')
    """
}
