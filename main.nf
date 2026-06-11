#!/usr/bin/env nextflow

// Imports
include { BOWTIE2_BUILD      } from './modules/nf-core/bowtie2/build/main'
include { ALIGN_PSEUDO_READS } from './modules/local/bowtie2/align/main'
include { MERGE_BAM          } from './modules/local/merge_bams/main'

// Define absolute paths to pseudo-reads and annotation
workflow {

    // Assemblies are in a directory
    assemblies_from_directory = channel.empty()
    if ( params.genomes_directory ){
        channel.fromPath("${params.genomes_directory}/*.{fna,fa,fasta}{,.gz}", checkIfExists: true)
        .map { f ->tuple(f.baseName.replaceFirst(/(\.fna|\.fa|\.fasta)(\.gz)?$/, ''),f)}
        .set { assemblies_from_directory }
    }

    // Assembly from a path
    assembly_from_path = channel.empty()
    if ( params.genome ){
        channel.fromPath(params.genome, checkIfExists: true)
        .map { f ->tuple(f.baseName.replaceFirst(/(\.fna|\.fa|\.fasta)(\.gz)?$/, ''),f)}
        .set { assembly_from_path }
    }

    // Mix
    assemblies_from_directory
        .mix( assembly_from_path )
        .unique()
        .set { assemblies_to_mask }

    // Pseudo-reads from directory
    if ( params.pseudo_reads_directory ){
        channel.fromPath("${params.pseudo_reads_directory}/*.{fna,fa,fasta}{,.gz}", checkIfExists: true)
        .set { pseudo_reads }
    }

    // Contig to species name file
    fna2name = channel.empty()
    if ( params.fna2name ){
        channel.fromPath(params.fna2name, checkIfExists: true)
        .set { fna2name }
    }

    // Create bowtie2 index for each assembly
    BOWTIE2_BUILD (
        assemblies_to_mask
    )

    // Prepare alignment input channel
    BOWTIE2_BUILD.out.index
        .combine( pseudo_reads )
        .map { id, index, reads ->
            [ id, index, reads, params.type_of_pseudo_reads, params.n_allowed_multimappers ]
        }
        .set { input_for_align }

    // Run alignment
    ALIGN_PSEUDO_READS (
        input_for_align
    )

    // Merge BAMs
    MERGE_BAM (
        ALIGN_PSEUDO_READS.out
        .groupTuple()
    )

    // Combine with fna2name
    MERGE_BAM.out.map { meta, bam, index ->
        [ meta, bam, index, params.type_of_pseudo_reads ]
    }
    .combine( fna2name )
    .set { detect_input }


    DETECT_EXOGENOUS(detect_input)

    make_bedfile(DETECT_EXOGENOUS.out.for_bedfile)

    make_bedfile.out.combine( assemblies_to_mask, by:0 ).view()

}




// Process 3: Detection
process DETECT_EXOGENOUS {

  publishDir params.outdir, mode: "copy"

  container 'docker://quay.io/biocontainers/mulled-v2-0697a5880de9863c66cba89c8310687052a940fc:c72ea422cf70582757ae5648f79b19857320259b-0'

  input:
    tuple val(input_ref), path(bam), path(bai), val(type_of_pseudo_reads), path(fna2name)

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
      ${type_of_pseudo_reads}

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
	[[ "\$f" == *.bam ]] && continue
    mv "\$f" "${input_ref}_${type_of_pseudo_reads}_\$f"
  done

  """
}

process make_bedfile {

  publishDir params.outdir, mode: "copy"

  input: 
    tuple val(ID), path(raw_bed)

  output: 
    tuple val(ID), path("*.bed")

  script:
  """
  out=\$(basename "$raw_bed" .txt).bed
  cut -f 2,3,4 "$raw_bed" | tail -n +2 | awk '{
    \$2=sprintf("%.0f",\$2);
    \$3=sprintf("%.0f",\$3);
    print}' OFS='\t' > "\$out"
  """
}

// mask the fasta with bedfile
process mask_fasta {

  conda 'bioconda::bedtools'

  publishDir params.outdir, mode: "copy"
  
  input: 
  tuple val(ID), path(bed), path(ref)

  output: path("*.masked.fna")

  script:
  """
    bedtools maskfasta -fi ${ref} -bed ${bed} -fo ${ID}.masked.fna
  """
}
