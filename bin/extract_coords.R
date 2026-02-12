#!/usr/bin/env Rscript

# FOR EXAMPLE: Rscript extract_coords.R GTDB
#
# AUTHOR: NIKOLAY OSKOLKOV, NBIS SCILIFELAB, EMAIL: nikolay.oskolkov@scilifelab.se

args = commandArgs(trailingOnly=TRUE)
TYPE_OF_PSEUDO_READS<-as.character(args[1])
myname<-as.character(args[2])#breadth of coverage file
fna2name<-as.character(args[3])

options(warn=0) #warnings are printed in error

contig<-strsplit(myname,"__")[[1]][1]
ref<-gsub(".boc","",strsplit(myname,"__")[[1]][2])
v<-as.numeric(readLines(myname))
r <- rle(v>0)
start <- cumsum(r$lengths)[r$values] - r$lengths[r$values] + 1
end <- start + r$lengths[r$values] - 1
coords<-data.frame(START=start,END=end)

output_contig<-data.frame(ORGANISM=rep(ref,dim(coords)[1]),CONTIG=rep(contig,dim(coords)[1]),START=coords$START,END=coords$END,LENGTH=coords$END-coords$START)


if(TYPE_OF_PSEUDO_READS=="RefSeq"){

average_depth<-vector(); n_reads<-vector(); micr1<-vector(); micr2<-vector(); micr3<-vector(); micr4<-vector(); micr5<-vector();
for(i in 1:dim(output_contig)[1])
{
average_depth<-append(average_depth,round(mean(v[output_contig$START[i]:output_contig$END[i]],na.rm=TRUE),0))
n_reads<-append(n_reads,system(paste0("samtools view PseudoReads_aligned_to_*.bam '",output_contig$CONTIG[i],"':",output_contig$START[i],"-",output_contig$END[i]," | wc -l"),intern=TRUE))
system("export LC_NUMERIC=en_US.utf-8")
micr<-system(paste0("samtools view PseudoReads_aligned_to_*.bam '",output_contig$CONTIG[i],"':",output_contig$START[i],"-",output_contig$END[i]," | cut -f1 | cut -f1,2 -d '_' | sort | uniq -c | sort -nr -k1,1 | head -5 | awk '{print $1\"_reads_\"$2}'"),intern=TRUE)
micr1<-append(micr1,micr[1]); micr2<-append(micr2,micr[2]); micr3<-append(micr3,micr[3]); micr4<-append(micr4,micr[4]); micr5<-append(micr5,micr[5])
}
n_reads<-as.numeric(n_reads)
output_contig$N_READS<-n_reads
output_contig$AVERAGE_DEPTH<-average_depth
output_contig$MICR1<-micr1; output_contig$MICR2<-micr2; output_contig$MICR3<-micr3; output_contig$MICR4<-micr4; output_contig$MICR5<-micr5
output_contig[is.na(output_contig)]<-"NA_reads_NA"

write.table(output_contig,file=paste0("coords_micr_like_regions_",ref,".txt"),col.names=!file.exists(paste0("coords_micr_like_regions_",ref,".txt")),append=T,row.names=FALSE,quote=FALSE,sep="\t")

}else if(TYPE_OF_PSEUDO_READS=="human"){

average_depth<-vector(); n_reads<-vector();
for(i in 1:dim(output_contig)[1])
{
average_depth<-append(average_depth,round(mean(v[output_contig$START[i]:output_contig$END[i]],na.rm=TRUE),0))
n_reads<-append(n_reads,system(paste0("samtools view PseudoReads_aligned_to_*.bam '",output_contig$CONTIG[i],"':",output_contig$START[i],"-",output_contig$END[i]," | wc -l"),intern=TRUE))
}
n_reads<-as.numeric(n_reads)
output_contig$N_READS<-n_reads
output_contig$AVERAGE_DEPTH<-average_depth

write.table(output_contig,file=paste0("coords_human_like_regions_",ref,".txt"),col.names=!file.exists(paste0("coords_human_like_regions_",ref,".txt")),append=T,row.names=FALSE,quote=FALSE,sep="\t")

}else{

average_depth<-vector(); n_reads<-vector(); micr1<-vector(); micr2<-vector(); micr3<-vector(); micr4<-vector(); micr5<-vector();
for(i in 1:dim(output_contig)[1])
{
average_depth<-append(average_depth,round(mean(v[output_contig$START[i]:output_contig$END[i]],na.rm=TRUE),0))
n_reads<-append(n_reads,system(paste0("samtools view PseudoReads_aligned_to_*.bam '",output_contig$CONTIG[i],"':",output_contig$START[i],"-",output_contig$END[i]," | wc -l"),intern=TRUE))
#system("export LC_NUMERIC=en_US.utf-8")
micr<-system(paste0("samtools view PseudoReads_aligned_to_*.bam '",output_contig$CONTIG[i],"':",output_contig$START[i],"-",output_contig$END[i]," | cut -f1 | sed 's/fna_/fna%/g' | cut -f1 -d '%' | sort | uniq -c | sort -nr -k1,1 | head -5 | awk '{print $1\"_reads_\"$2}'"),intern=TRUE)
micr1<-append(micr1,micr[1]); micr2<-append(micr2,micr[2]); micr3<-append(micr3,micr[3]); micr4<-append(micr4,micr[4]); micr5<-append(micr5,micr[5])
}
n_reads<-as.numeric(n_reads)
output_contig$N_READS<-n_reads
output_contig$AVERAGE_DEPTH<-average_depth
output_contig$MICR1<-micr1; output_contig$MICR2<-micr2; output_contig$MICR3<-micr3; output_contig$MICR4<-micr4; output_contig$MICR5<-micr5
output_contig[is.na(output_contig)]<-"NA_reads_NA"

annot<-read.delim(fna2name,header=FALSE,sep="\t")

micr1_sep<-as.data.frame(matrix(unlist(strsplit(as.character(output_contig$MICR1),"_reads_")),ncol=2,byrow=TRUE))
micr1_sep$V3<-annot$V3[match(as.character(micr1_sep$V2),as.character(annot$V1))]
micr1_sep$V3<-gsub(" ","_",as.character(micr1_sep$V3))
output_contig$MICR1<-paste0(micr1_sep$V1,"_reads_",micr1_sep$V3)

micr2_sep<-as.data.frame(matrix(unlist(strsplit(as.character(output_contig$MICR2),"_reads_")),ncol=2,byrow=TRUE))
micr2_sep$V3<-annot$V3[match(as.character(micr2_sep$V2),as.character(annot$V1))]
micr2_sep$V3<-gsub(" ","_",as.character(micr2_sep$V3))
output_contig$MICR2<-paste0(micr2_sep$V1,"_reads_",micr2_sep$V3)

micr3_sep<-as.data.frame(matrix(unlist(strsplit(as.character(output_contig$MICR3),"_reads_")),ncol=2,byrow=TRUE))
micr3_sep$V3<-annot$V3[match(as.character(micr3_sep$V2),as.character(annot$V1))]
micr3_sep$V3<-gsub(" ","_",as.character(micr3_sep$V3))
output_contig$MICR3<-paste0(micr3_sep$V1,"_reads_",micr3_sep$V3)

micr4_sep<-as.data.frame(matrix(unlist(strsplit(as.character(output_contig$MICR4),"_reads_")),ncol=2,byrow=TRUE))
micr4_sep$V3<-annot$V3[match(as.character(micr4_sep$V2),as.character(annot$V1))]
micr4_sep$V3<-gsub(" ","_",as.character(micr4_sep$V3))
output_contig$MICR4<-paste0(micr4_sep$V1,"_reads_",micr4_sep$V3)

micr5_sep<-as.data.frame(matrix(unlist(strsplit(as.character(output_contig$MICR5),"_reads_")),ncol=2,byrow=TRUE))
micr5_sep$V3<-annot$V3[match(as.character(micr5_sep$V2),as.character(annot$V1))]
micr5_sep$V3<-gsub(" ","_",as.character(micr5_sep$V3))
output_contig$MICR5<-paste0(micr5_sep$V1,"_reads_",micr5_sep$V3)

write.table(output_contig,file=paste0("coords_micr_like_regions_",ref,".txt"),col.names=!file.exists(paste0("coords_micr_like_regions_",ref,".txt")),append=T,row.names=FALSE,quote=FALSE,sep="\t")

}

