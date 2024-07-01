#!/bin/bash


# Usage: ./RNA_flow_alignment_directRNA_ONT.sh sample_name

# Tasks:
#   	1. convert RNA --> DNA
#       2. align with minimap2
#       3. filter uniquely mapped reads
#
#
##########################################################################################
##########################################################################################
##########################################################################################
#

module load gcc/6.2.0
module load samtools/1.3.1

sample_name=$1


###############################################

#       1. Convert RNA sequence to DNA

mkdir -p /path/to/startingFiles/$sample_name/RNAtoDNA
outDir=/path/to/startingFiles/$sample_name/RNAtoDNA
inDir=/path/to/startingFiles/$sample_name/fastq_pass

# Combine all fastq files into one
for i in $inDir/*.fastq.gz ; do zcat $i >> ${inDir}/${sample_name}_pass.fastq ; done

# Convert Us to Ts
awk '{if(NR%4==2){gsub(/U/,"T",$0); print} else {print $0}}' ${inDir}/${sample_name}_pass.fastq > ${outDir}/${sample_name}_pass_RNAtoDNA.fastq


#############################################

#       2. Run minimap2 alignment (with ONT default parameters) for all samples


programDir=/path/to/software/minimap2
reference=/path/to/genomes/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.and.sacCer3.fa

mkdir -p /path/to/alignment/${sample_name}/minimap2/LogErr
outDir=/path/to/alignment/${sample_name}/minimap2
inDir=/path/to/startingFiles/${sample_name}/RNAtoDNA
reads=${inDir}/${sample_name}_pass_RNAtoDNA.fastq


${programDir}/minimap2 -ax splice -uf -k14 ${reference} ${reads} > ${outDir}/${sample_name}_minimap2.sam
samtools sort ${outDir}/${sample_name}_minimap2.sam -o ${outDir}/${sample_name}_minimap2_sort.bam
samtools index ${outDir}/${sample_name}_minimap2_sort.bam


###############################################3


# 3. Extract unique reads

Dir=/path/to/alignment/${sample_name}/minimap2

cd ${Dir}

grep ^@ ${sample_name}_minimap2.sam > headers.sam
grep -v ^@  ${sample_name}_minimap2.sam | sort > alignment.sam
awk '$3!="*" {print}' alignment.sam | cut -f 1 | sort | uniq -c | awk '$1==1 {print $2}' > uniq_names.txt
grep -F -f uniq_names.txt alignment.sam > uniq_alignment.sam
cat headers.sam uniq_alignment.sam > headers_uniq_temp.sam
samtools view -bT ${reference} headers_uniq_temp.sam > headers_uniq_temp.bam
samtools sort headers_uniq_temp.bam -o ${sample_name}_minimap2_uniq_sort.bam
samtools index ${sample_name}_minimap2_uniq_sort.bam
rm headers.sam
rm alignment.sam
rm uniq_names.txt
rm uniq_alignment.sam
rm headers_uniq_temp.sam
rm headers_uniq_temp.bam

