#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define absolute paths to pseudo-reads and annotation
workflow {

    pseudo_reads_file = Channel.fromPath("${params.pseudo_reads_file_dir}/*.{fna,fa,fasta}{,.gz}")
    .ifEmpty {
        log.error "No input files found in ${params.pseudo_reads_file_dir}"
        System.exit(1)
    }

    output_dir = Channel.fromPath(params.output_dir)

    files = Channel.fromPath("${params.input_dir}/*.{fna,fa,fasta}{,.gz}")
    .ifEmpty {
        log.error "No input files (.fna, .fa, .fasta, optionally .gz) found in ${params.input_dir}"
        System.exit(1)
    }
    .map { f ->
        tuple(
            f.baseName.replaceFirst(/(\.fna|\.fa|\.fasta)(\.gz)?$/, ''),
            f
            )
    }
    .view()

    index_reference(files)
    
    input_for_align = index_reference.out
     			.combine(pseudo_reads_file)
     			.map{it -> tuple(it, params.type_of_pseudo_reads, params.n_allowed_multimappers)}
                .map{it -> it.flatten()}
                .map{it -> tuple(it[0], tuple(it[1],it[2],it[3],it[4],it[5],it[6]),it[7], params.type_of_pseudo_reads, params.n_allowed_multimappers)}

    input_for_align.view()

    align_pseudo_reads(input_for_align)

    merge_bam(align_pseudo_reads.out.groupTuple())

    merge_bam.out
    .map{it -> tuple(it, params.type_of_pseudo_reads, params.work_dir, params.fna2name)}
    .map{it -> it.flatten()}
    .set{  detect_input }

    // detect_exogenous(align_pseudo_reads.out.bam, align_pseudo_reads.out.ref, params.type_of_pseudo_reads, params.work_dir, params.fna2name)
    detect_exogenous(detect_input)

    make_bedfile(detect_exogenous.out.for_bedfile)

    make_bedfile.out.combine( files, by:0 ).view()
    // mask_fasta(make_bedfile.out.combine( files,by:0 ))

}

// Process 1: Indexing
process index_reference {
    conda './envs/bowtie2.yml'

    cpus { 40 * task.attempt }
    memory { 32.GB * task.attempt }
    time { 1.hour * task.attempt }
        
    errorStrategy {
    if( task.exitStatus == null || task.exitStatus in 137..140 || task.exitStatus == 143 ) {
        // - represent external termination due to time limit
        return task.attempt <= 6 ? 'retry' : 'terminate'
    }
    else {
        return 'terminate'
    }
    }


    input:
    tuple val(ID), path(input_ref)

    output:
    tuple val(ID), path("*.bt2l")

    script:
    """
    bowtie2-build --large-index \$(basename ${input_ref}) \$(basename ${input_ref}) --threads "${task.cpus}"
    """
}

// Process 2: Alignment
process align_pseudo_reads {
    conda './envs/bowtie2.yml'

    cpus { 10 * task.attempt }
    memory { 8.GB * task.attempt }
    time { 30.m * task.attempt }

    errorStrategy {
    if( task.exitStatus == null || task.exitStatus in 137..140 || task.exitStatus == 143 ) {
        return task.attempt <= 6 ? 'retry' : 'terminate'
    }
    else {
        return 'terminate'
    }
    }

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


// merge bam of same sample mapped to different databases
process merge_bam {

    cpus { 25 * task.attempt }
    memory { 20.GB * task.attempt }
    time { 1.hour * task.attempt }

    errorStrategy {
    if( task.exitStatus == null || task.exitStatus in 137..140 || task.exitStatus == 143 ) {
        return task.attempt <= 6 ? 'retry' : 'terminate'
    }
    else {
        return 'terminate'
    }
    }
    
    conda './envs/bowtie2.yml'
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


// Process 3: Detection
process detect_exogenous {

  publishDir params.output_dir, 
        mode: "copy"

  cpus { 20 * task.attempt }
  memory { 15.GB * task.attempt }
  time { 4.hour * task.attempt }

    errorStrategy {
    if( task.exitStatus == null || task.exitStatus in 137..140 || task.exitStatus == 143 ) {
        return task.attempt <= 6 ? 'retry' : 'terminate'
    }
    else {
        return 'terminate'
    }
    }

  container 'docker://quay.io/biocontainers/mulled-v2-0697a5880de9863c66cba89c8310687052a940fc:c72ea422cf70582757ae5648f79b19857320259b-0'

  input:
    tuple val(input_ref), path(bam), path(bai), val(type_of_pseudo_reads), val(work_dir), path(fna2name)

  output:
    path("*abund_*.txt")
    tuple val(input_ref), path("*coords_micr_like_regions*.txt"), emit: for_bedfile
    path("*boc_*.txt")
    path("*_microbes_abundant_*.txt")
    
    

  script:
  """
  #get just bam
  bamfile=\$(echo $bam | awk '{print \$1}' )

  detect_exogenous.sh \
      \${bamfile} \
      ${input_ref} \
      ${type_of_pseudo_reads} \
      ${work_dir}

  echo "GENERATE COORDINATIONS OF MICROBIAL-LIKE REGIONS (BEDFILES)"
  for j in \$(cat refs_uniq_sorted.txt)
	do
	echo \${j} CONTIG OF ${input_ref}
	extract_coords.R ${type_of_pseudo_reads} \${j}__${input_ref}.boc $fna2name
	echo DELETING BAM AND COMPRESSING BOC FILES
	rm \${j}.bam
	rm \${j}__${input_ref}.boc
  done
  #remove intermediate files
  rm refs_uniq_sorted.txt refs_uniq_sorted_reads.txt total_length_per_ref.txt boc_per_ref.txt
  rm $fna2name #avoid output it

  #add prefix
  for f in \$(ls * | grep -v .bam); do
    mv "\$f" "${input_ref}_${type_of_pseudo_reads}_\$f"
  done

  """
}

process make_bedfile {

  publishDir params.output_dir, mode: "copy"

  cpus { 1 * task.attempt }
  memory { 0.8.GB * task.attempt }
  time { 20.m * task.attempt }

    errorStrategy {
    if( task.exitStatus == null || task.exitStatus in 137..140 || task.exitStatus == 143 ) {
        return task.attempt <= 6 ? 'retry' : 'terminate'
    }
    else {
        return 'terminate'
    }
    }

  input: 
    tuple val(ID), path(raw_bed)

  output: 
    tuple val(ID), path("*.bed")

  script:
  """
  out=\$(basename "$raw_bed" .txt).bed
  cut -f 2,3,4 "$raw_bed" | tail -n +2 > "\$out"
  """
}

// mask the fasta with bedfile
process mask_fasta {

  conda 'bioconda::bedtools'

  publishDir params.output_dir, mode: "copy"

  cpus { 1 * task.attempt }
  memory { 0.8.GB * task.attempt }
  time { 20.m * task.attempt }

    errorStrategy {
    if( task.exitStatus == null || task.exitStatus in 137..140 || task.exitStatus == 143 ) {
        return task.attempt <= 6 ? 'retry' : 'terminate'
    }
    else {
        return 'terminate'
    }
    }

  input: 
  tuple val(ID), path(bed), path(ref)

  output: path("*.masked.fna")

  script:
  """
    bedtools maskfasta -fi ${ref} -bed ${bed} -fo ${ID}.masked.fna
  """
}