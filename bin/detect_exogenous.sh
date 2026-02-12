#!/bin/bash

# AUTHOR: NIKOLAY OSKOLKOV, NBIS SCILIFELAB, EMAIL: nikolay.oskolkov@scilifelab.se

BAM=$1
REF_GENOME=$2
TYPE_OF_PSEUDO_READS=$3
WORK_DIR=$4

printf "\n"; echo RANKING ${REF_GENOME} CONTIGS BY NUMBER OF MAPPED MICROBIAL READS
samtools view $BAM | cut -f3 | sort | uniq -c | sort -nr -k1,1 | awk '{ t = $1; $1 = $2; $2 = t; print; }' | tr ' ' '\t' > contigs_abund_sorted_${TYPE_OF_PSEUDO_READS}_${REF_GENOME}.txt
awk '{print$1}' contigs_abund_sorted_${TYPE_OF_PSEUDO_READS}_${REF_GENOME}.txt > refs_uniq_sorted.txt
awk '{print$2}' contigs_abund_sorted_${TYPE_OF_PSEUDO_READS}_${REF_GENOME}.txt > refs_uniq_sorted_reads.txt
sed -i '1s/^/CONTIG\tN_READS\n/' contigs_abund_sorted_${TYPE_OF_PSEUDO_READS}_${REF_GENOME}.txt
sed -i "1s/^/ORGANISM\t/; 2,\$s/^/${REF_GENOME}\t/" contigs_abund_sorted_${TYPE_OF_PSEUDO_READS}_${REF_GENOME}.txt

echo COMPUTING BREADTH OF COVERAGE FOR EACH CONTIG AND COORDINATES OF MICROBIAL CONTAMINATION FOR ${REF_GENOME} REFERENCE GENOME
for j in $(cat refs_uniq_sorted.txt)
do
echo ${j} CONTIG OF ${REF_GENOME}
samtools view -b $BAM ${j} > ${j}.bam
samtools depth -g 0x100 -a ${j}.bam | cut -f3 > ${j}__${REF_GENOME}.boc
awk -v covered_length=$(awk '{if($1>0)print$0}' ${j}__${REF_GENOME}.boc | wc -l) -v total_length=$(wc -l ${j}__${REF_GENOME}.boc | cut -f1 -d ' ') 'BEGIN { print ( covered_length / total_length ) }' >> boc_per_ref.txt
wc -l ${j}__${REF_GENOME}.boc | cut -f1 -d ' ' >> total_length_per_ref.txt
#echo EXTRACTING COORDINATES OF MICROBIAL CONTAMINATION
#Rscript $WORK_DIR/bin/extract_coords.R ${TYPE_OF_PSEUDO_READS}
#echo DELETING BAM AND COMPRESSING BOC FILES
#rm ${j}.bam
#gzip ${j}__${REF_GENOME}.boc
done

printf "\n"; echo AGGREGATING RESULTS FOR ${REF_GENOME} REFERENCE GENOME AND CLEANING
#export LC_NUMERIC=en_US.utf-8
paste refs_uniq_sorted.txt refs_uniq_sorted_reads.txt total_length_per_ref.txt boc_per_ref.txt | sort -gr -k4,4 > contigs_boc_sorted_${TYPE_OF_PSEUDO_READS}_${REF_GENOME}.txt
sed -i '1s/^/CONTIG\tN_READS\tLENGTH\tBOC\n/' contigs_boc_sorted_${TYPE_OF_PSEUDO_READS}_${REF_GENOME}.txt
sed -i "1s/^/ORGANISM\t/; 2,\$s/^/${REF_GENOME}\t/" contigs_boc_sorted_${TYPE_OF_PSEUDO_READS}_${REF_GENOME}.txt
#rm refs_uniq_sorted.txt refs_uniq_sorted_reads.txt total_length_per_ref.txt boc_per_ref.txt

echo COMPUTING LIST OF MOST ABUNDANT MICROBES CONTMAMINATING ${REF_GENOME} REFERENCE GENOME
if [ ${TYPE_OF_PSEUDO_READS} = RefSeq ]; then
    samtools view $BAM | cut -f1 | cut -f1,2 -d '_' | sort | uniq -c | sort -nr -k1,1 > microbes_abundant_${TYPE_OF_PSEUDO_READS}_${REF_GENOME}.txt
elif [ ${TYPE_OF_PSEUDO_READS} = human ]; then
    echo "NO LIST OF ABUNDANT ORGANISMS GENERATED BECAUSE ALL PSEUDO-READS ARE HUMAN"
else
    samtools view $BAM | cut -f1 | sed 's/fna_/fna%/g' | cut -f1 -d '%' | sort | uniq -c | sort -nr -k1,1 > microbes_abundant_${TYPE_OF_PSEUDO_READS}_${REF_GENOME}.txt
fi

printf "\n"; echo ANALYSIS FOR ${REF_GENOME} REFERENCE GENOME FINISHED SUCCESSFULLY
