# THIS WORKFLOW DETECTS REGIONS OF POTENTIAL MICROBIAL CONTAMINATION IN PROVIDED EUKARYOTIC REFERENCE GENOME. IT ACCEPTS A REFERENCE GENOME IN FASTA-FORMAT AND OUTPUTS A BED-FILE WITH COORDINATES OF MICROBIAL CONTAMINATION.
#
# THE WORKFLOW REQUIRES THE HELPING R-SCRIPT extract_coords_micr_contam.R TO BE PLACED NEXT TO micr_cont_detect.sh
# THE WORKFLOW REQUIRES R, AWK, BOWTIE2 AND SAMTOOLS TO BE INSTALLED AND AVAILABE IN THE PATH
#
# RUN THE WORKFLOW AS:  ./micr_cont_detect.sh REF_GENOME INPUT_DIR REFSEQ_OR_GTDB THREADS MICR_READS GTDB_ANNOT
# WHERE: REF_GENOME     - GZIPPED FASTA-REFERENCE
#        INPUT_DIR      - DIRECTORY CONTAINING THE FASTA-REFERENCE
#        REFSEQ_OR_GTDB - WHETHER REFSEQ OR GTDB SLICED MICROBIAL PSEUDO-READS SHOULD BE USED FOR DETECTION OF CONTAMINATED REGIONS IN THE FASTA-REFERENCE, CAN BE "RefSeq" OR "GTDB"
#        THERADS        - NUMBER OF THERADS AVAILABLE
#        MICR_READS     - GTDB OR REFSEQ MICROBIAL PSEUDO-READS (GZIPPED FASTA-FILE) PROVIDED TOGETHER WITH THE WORKFLOW
#        GTDB_ANNOT     - GTDB ANNOTATION FILE GTDB_fna2name.txt PROVIDED TOGETHER WITH THE WORKFLOW
#
# FOR EXAMPLE: ./micr_cont_detect.sh GCA_027887165.1_mMonDom1.pri_genomic.fna.gz /home/data GTDB 20 GTDB_sliced_seqs_sliding_window.fna.gz GTDB_fna2name.txt
#
# AUTHOR: NIKOLAY OSKOLKOV, NBIS SCILIFELAB, EMAIL: nikolay.oskolkov@scilifelab.se

#!/bin/bash

REF_GENOME=$1
INPUT_DIR=$2
REFSEQ_OR_GTDB=$3
THREADS=$4
MICR_READS=$5
GTDB_ANNOT=$6


printf "\n"; echo PREPARING FILES FOR ANALYSIS OF ${REF_GENOME} REFERENCE GENOME
cd "${INPUT_DIR}"
mkdir ${REF_GENOME}_${REFSEQ_OR_GTDB}

printf "\n"; echo BUILDING BOWTIE2 INDEX FOR ${REF_GENOME} REFERENCE GENOME
bowtie2-build --large-index ${REF_GENOME} ${REF_GENOME} --threads ${THREADS} >> bowtie2-build.log 2>&1
cd ${REF_GENOME}_${REFSEQ_OR_GTDB}

echo ALIGNING MICROBIAL READS WITH BOWTIE2 TO ${REF_GENOME} REFERENCE GENOME
bowtie2 --large-index -f -k 10 -x ../${REF_GENOME} --end-to-end --quiet --threads ${THREADS} --very-sensitive -U ../../${MICR_READS} | samtools view -bS -F 4 -h -@ ${THREADS} - | samtools sort -@ ${THREADS} - > MicrReads_aligned_to_${REF_GENOME}.bam
samtools index -c MicrReads_aligned_to_${REF_GENOME}.bam

printf "\n"; echo RANKING ${REF_GENOME} CONTIGS BY NUMBER OF MAPPED MICROBIAL READS
samtools view MicrReads_aligned_to_${REF_GENOME}.bam | cut -f3 | sort | uniq -c | sort -nr -k1,1 | awk '{ t = $1; $1 = $2; $2 = t; print; }' | tr ' ' '\t' > contigs_abund_sorted_${REFSEQ_OR_GTDB}_${REF_GENOME}.txt
awk '{print$1}' contigs_abund_sorted_${REFSEQ_OR_GTDB}_${REF_GENOME}.txt > refs_uniq_sorted.txt
awk '{print$2}' contigs_abund_sorted_${REFSEQ_OR_GTDB}_${REF_GENOME}.txt > refs_uniq_sorted_reads.txt
sed -i '1s/^/CONTIG\tN_READS\n/' contigs_abund_sorted_${REFSEQ_OR_GTDB}_${REF_GENOME}.txt
sed -i "1s/^/ORGANISM\t/; 2,\$s/^/${REF_GENOME}\t/" contigs_abund_sorted_${REFSEQ_OR_GTDB}_${REF_GENOME}.txt

echo COMPUTING BREADTH OF COVERAGE FOR EACH CONTIG AND COORDINATES OF MICROBIAL CONTAMINATION FOR ${REF_GENOME} REFERENCE GENOME
for j in $(cat refs_uniq_sorted.txt)
do
echo ${j} CONTIG OF ${REF_GENOME}
samtools view -b MicrReads_aligned_to_${REF_GENOME}.bam ${j} > ${j}.bam
samtools depth -g 0x100 -a ${j}.bam | cut -f3 > ${j}__${REF_GENOME}.boc
awk -v covered_length=$(awk '{if($1>0)print$0}' ${j}__${REF_GENOME}.boc | wc -l) -v total_length=$(wc -l ${j}__${REF_GENOME}.boc | cut -f1 -d ' ') 'BEGIN { print ( covered_length / total_length ) }' >> boc_per_ref.txt
wc -l ${j}__${REF_GENOME}.boc | cut -f1 -d ' ' >> total_length_per_ref.txt
echo EXTRACTING COORDINATES OF MICROBIAL CONTAMINATION
Rscript ../../extract_coords_micr_contam.R ${REFSEQ_OR_GTDB} ../../${GTDB_ANNOT}
echo DELETING BAM AND COMPRESSING BOC FILES
rm ${j}.bam
gzip ${j}__${REF_GENOME}.boc
done

printf "\n"; echo AGGREGATING RESULTS FOR ${REF_GENOME} REFERENCE GENOME AND CLEANING
export LC_NUMERIC=en_US.utf-8
paste refs_uniq_sorted.txt refs_uniq_sorted_reads.txt total_length_per_ref.txt boc_per_ref.txt | sort -gr -k4,4 > contigs_boc_sorted_${REFSEQ_OR_GTDB}_${REF_GENOME}.txt
sed -i '1s/^/CONTIG\tN_READS\tLENGTH\tBOC\n/' contigs_boc_sorted_${REFSEQ_OR_GTDB}_${REF_GENOME}.txt
sed -i "1s/^/ORGANISM\t/; 2,\$s/^/${REF_GENOME}\t/" contigs_boc_sorted_${REFSEQ_OR_GTDB}_${REF_GENOME}.txt
rm refs_uniq_sorted.txt refs_uniq_sorted_reads.txt total_length_per_ref.txt boc_per_ref.txt

echo COMPUTING LIST OF MOST ABUNDANT MICROBES CONTMAMINATING ${REF_GENOME} REFERENCE GENOME
if [ ${REFSEQ_OR_GTDB} = RefSeq ]; then
    samtools view MicrReads_aligned_to_${REF_GENOME}.bam | cut -f1 | cut -f1,2 -d '_' | sort | uniq -c | sort -nr -k1,1 > microbes_abundant_${REFSEQ_OR_GTDB}_${REF_GENOME}.txt
else
    samtools view MicrReads_aligned_to_${REF_GENOME}.bam | cut -f1 | sed 's/fna_/fna%/g' | cut -f1 -d '%' | sort | uniq -c | sort -nr -k1,1 > microbes_abundant_${REFSEQ_OR_GTDB}_${REF_GENOME}.txt
fi

printf "\n"; echo ANALYSIS FOR ${REF_GENOME} REFERENCE GENOME FINISHED SUCCESSFULLY
