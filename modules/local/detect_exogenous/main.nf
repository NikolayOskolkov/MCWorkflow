process DETECT_EXOGENOUS {

  label 'process_low'
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