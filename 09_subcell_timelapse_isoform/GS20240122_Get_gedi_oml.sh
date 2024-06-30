#!/bin/bash
#Robert Ietswaart 
#source: https://github.com/erhard-lab/gedi/wiki/Preparing-genomes

date=20240126
baseDir=/n/groups/churchman/ri23/bseq/GS${date}

program=/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi 

scDir=/n/groups/churchman/ri23/code/bseq/

module load gcc/9.2.0
module load java/jdk-1.8u112
# module load bedops/2.4.30
module load bedtools/2.30.0
module load samtools/1.15.1
module load R/4.3.1

cd ${baseDir}
mkdir -p ${baseDir}/LogErr


#Source: GS20240122_Get_gene_level_gtf.sh 
gtf=K562_ensGRCh38_MTmod_dm6_ercc_cat_genes_only.gtf


fa=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_dm6_ercc_cat.fasta


#Generate a GRAND-SLAM index file (.oml) format
echo "begin 1"
$program -e IndexGenome -s $fa -a $gtf -nobowtie -nokallisto -nostar -p -D 
#unnecessary arguments
#-D output debugging info
#-p output progress

#This generates the .oml file: 
#/home/ri23/.gedi/genomic/K562_ensGRCh38_MTmod_dm6_ercc_cat_genes_only.oml
cp /home/ri23/.gedi/genomic/K562_ensGRCh38_MTmod_dm6_ercc_cat_genes_only.oml ./

echo "end of script"