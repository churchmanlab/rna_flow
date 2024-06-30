#!/bin/bash
#Robert Ietswaart 
#Get introns_sort.bed which is necessary for KDii sample processing in 


date=20240126
baseDir=/n/groups/churchman/ri23/bseq/GS${date}/
scDir=/n/groups/churchman/ri23/code/bseq/

module load gcc/9.2.0
module load java/jdk-11.0.11 #jdk-1.8u112
# module load bedops/2.4.30
module load bedtools/2.30.0
# module load samtools/1.15.1
module load R/4.3.1

cd ${baseDir}
mkdir -p ${baseDir}LogErr

fai=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_dm6_ercc_cat.fasta.fai
gtf=K562_ensGRCh38_MTmod_dm6_ercc_cat.gtf  #From GS20240122_Get_gene_level_gtf.sh

#1. Sort gtf file
sort -k1,1 -k4,4n -k5,5n $gtf > K562_ensGRCh38_MTmod_dm6_ercc_cat_sort.gtf

#2. generate genome.file that contains the chromosome start and end site
awk '{print $1"\t"$2}' $fai > genome.file
sort -k1,1 -k2,2n genome.file > genome_sort.file

#3 Filter for only exon features
awk 'BEGIN{FS=OFS="\t"} $3 == "exon"' K562_ensGRCh38_MTmod_dm6_ercc_cat_sort.gtf > exons.gtf
sort -k1,1 -k4,4n -k5,5n exons.gtf > exons_sort.gtf

#4 Generate a bed file of purely intronic regions, i.e. there is no transcript isoform that contains
bedtools complement -i exons.gtf -g genome_sort.file > introns.bed
sort -k1,1 -k2,2n -k3,3n introns.bed > introns_sort.bed

#5. Get intron.gtf
Rscript ${scDir}GS20240122_Get_introns_gtf.R exons.gtf ${baseDir}
sort -k1,1 -k4,4n -k5,5n introns.gtf > introns_sort.gtf


# #Help page
# bedtools complement -h