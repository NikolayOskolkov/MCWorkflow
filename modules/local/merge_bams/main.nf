
process MERGE_BAM {

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    tuple val(ID), path(bams)

    output:
    tuple val(ID), path("*_merged.sorted.bam"), path("*.bam.csi")

    script:
    """
        #filtering out unmapped reads in case it's not done for input bam
        for bam1 in *.bam; do
            [[ "\$bam1" == *mapped.bam ]] && continue
            samtools view -@ "${task.cpus}" -b -F 0x4 "\$bam1" -o "\$(basename \$bam1 .bam).mapped.bam"
        done
    
        samtools merge ${ID}.merged.bam *.mapped.bam

        samtools quickcheck ${ID}.merged.bam || {
            echo "ERROR: Merging is not successful: ${ID}.merged.bam" >&2
            exit 1
        }

        samtools sort -@ "${task.cpus}" -o ${ID}_merged.sorted.bam ${ID}.merged.bam

        samtools index -c ${ID}_merged.sorted.bam
    """
}
