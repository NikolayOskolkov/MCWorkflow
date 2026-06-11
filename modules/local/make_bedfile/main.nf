
process MAKE_BEDFILE {

    label 'process_low'

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
