#!/usr/bin/env nextflow

// Imports
include { BOWTIE2_BUILD      } from './modules/nf-core/bowtie2/build/main'
include { ALIGN_PSEUDO_READS } from './modules/local/bowtie2/align/main'
include { MERGE_BAM          } from './modules/local/merge_bams/main'
include { DETECT_EXOGENOUS   } from './modules/local/detect_exogenous/main'
include { MAKE_BEDFILE       } from './modules/local/make_bedfile/main'

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

    // Run the detect exogenous scripts
    DETECT_EXOGENOUS(
        detect_input
    )

    // Create bedfiles
    MAKE_BEDFILE(
        DETECT_EXOGENOUS.out.for_bedfile
    )

}
