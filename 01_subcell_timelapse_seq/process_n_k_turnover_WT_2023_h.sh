#!/bin/bash

#SBATCH -c 4
#SBATCH -t 60
#SBATCH -p short
#SBATCH --mem=25G
#SBATCH -o %j.out
#SBATCH -e %j.err

## To be run following initial grand slam output, before running findMismatches to generate [n,k] matrix
## Make sure that GS has been run a second time without specifying -no4sUpattern so the unlabeled library's MAPs are reported in the GS output

# SET sample name - can also set manually (sample, timepoint, experiment, experiment 2)
sample=`pwd | cut -d '/' -f 8`
timepoint=2
x=`echo $sample | cut -c 1 `
y=`echo $sample | cut -d '_' -f 2 `
experiment=${x}-${y}
experiment2=${x}_${y}

# SET path to raw GS output
mini_path=`pwd | cut -d '/' -f 6`
gs_out=/n/groups/churchman/bms36/${mini_path}/GS_2023/${experiment}/${experiment2}.tsv

# SET path to list of ptc genes for that genome
ptc_list=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_MTmod_ptc_list.txt
#ptc_list=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/NIH3T3_mm10_MTmod_ptc_list.txt
#OLD ptc_list=/n/groups/churchman/bms36/sandbox/TimelapseSeq/human_ptc_list.txt
#OLD ptc_list=/n/groups/churchman/bms36/sandbox/TimelapseSeq/mouse_ptc_list.txt

# SET path to genome gtf 
gtf_list=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_MTmod.noHeader.withIDs.gtf
#gtf_list=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/NIH3T3_mm10_MTmod.noHeader.withIDs.gtf
#OLD gtf_list=/n/groups/churchman/bms36/genomes/hg38_dm6_ercc_grandslam/mtc_snp_masked/K562_ensGRCh38_MTmod_dm6_ercc_cat_hg_only_no_MT.noHeader.withID.gtf
#OLD gtf_list=/n/groups/churchman/bms36/genomes/mm10_dm6_ercc_grandslam/mtc_snp_masked/mouseNIH3T3_mm10_MTmod_dm6_ercc_cat_mm10_only_no_MT.noHeader.withID.gtf

# SET path to genome fasta
genome_fasta=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_dm6_ercc_cat.fasta
#genome_fasta=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/NIH3T3_mm10_dm6_ercc_cat.fasta
#OLD genome_fasta=/n/groups/churchman/bms36/genomes/mm10_dm6_ercc_grandslam/mtc_snp_masked/mouseNIH3T3_mm10_dm6_ercc_cat.fasta
#OLD genome_fasta=/n/groups/churchman/bms36/genomes/hg38_dm6_ercc_grandslam/mtc_snp_masked/K562_ensGRCh38_dm6_ercc_cat.fasta

# load required modules
module load gcc/9.2.0
module load samtools/1.15.1
module load java/jdk-11.0.11
module load bedtools/2.30.0
module load R/4.3.1

# run top100genes_turnover.R
/n/groups/churchman/bms36/programs/top1000genes_turnover.R $gs_out $ptc_list $gtf_list $sample $timepoint

# get reads that align to the top 100 genes
samtools view -@ 3 -b -L ${sample}_top1000genes_turnover.bed ${sample}_sort.human_noMT.bam > ${sample}_top1000genes_turnover.bam

#echo "Running samtools flagstat: all input reads"
samtools flagstat -@ 3 ${sample}_sort.human_noMT.bam

#echo "Running samtools flagstat: top 100 genes"
samtools flagstat -@ 3 ${sample}_top1000genes_turnover.bam

# convert from .bam to .sam
samtools view -@ 3 -h -O SAM -o ${sample}_top1000genes.paired.sam ${sample}_top1000genes_turnover.bam
rm ${sample}_top1000genes_turnover.bam

# remove soft clipped bases from reads
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/biostar84452.jar -o ${sample}_top1000genes.paired.noSoft.sam ${sample}_top1000genes.paired.sam
rm ${sample}_top1000genes.paired.sam

# converting cigar string to discriminate between mismatches and matches
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/samfixcigar.jar -r $genome_fasta ${sample}_top1000genes.paired.noSoft.sam > ${sample}_top1000genes.paired.noSoft.fixCigar.sam
rm ${sample}_top1000genes.paired.noSoft.sam

# get rid of sam header
samtools view -@ 3 -O SAM -o ${sample}_top1000genes.paired.noSoft.fixCigar.noHead.sam ${sample}_top1000genes.paired.noSoft.fixCigar.sam

# convert to bam (the sam with header)
samtools view -@ 3 -b -o ${sample}_top1000genes.paired.noSoft.fixCigar.bam ${sample}_top1000genes.paired.noSoft.fixCigar.sam
rm ${sample}_top1000genes.paired.noSoft.fixCigar.sam

# convert to modified bed format
bedtools bamtobed -cigar -i ${sample}_top1000genes.paired.noSoft.fixCigar.bam > ${sample}_reads_top1000genes.bed
rm ${sample}_top1000genes.paired.noSoft.fixCigar.bam

# run modifyBed.R script to edit bed file so MM analysis can be run on it *Need more than 25G memory for this step
/n/groups/churchman/bms36/programs/modifyBed.R ${sample}_reads_top1000genes.bed ${sample}_top1000genes.paired.noSoft.fixCigar.noHead.sam
rm ${sample}_reads_top1000genes.bed
rm ${sample}_top1000genes.paired.noSoft.fixCigar.noHead.sam

# make sure the resulting file from modifyBed.R contains same read numbers as above
echo "Number of reads in top1000 bed file, before removing singles:"
echo $(wc -l ${sample}_reads_top1000genes.bedwithC.bed)

# keep only reads that are still paired correctly
cut -f 4 ${sample}_reads_top1000genes.bedwithC.bed | awk -F'/' 'BEGIN{OFS="\t";} {print $1}' > readPrefix.txt
paste ${sample}_reads_top1000genes.bedwithC.bed readPrefix.txt > readsPlusPrefix.txt
rm readPrefix.txt
rm ${sample}_reads_top1000genes.bedwithC.bed
/n/groups/churchman/bms36/sandbox/TimelapseSeq/MismatchScripts_batchSubmission/sortReads_byName_noSingles.R readsPlusPrefix.txt
rm readsPlusPrefix.txt
uniq -f 10 -D readsPlusPrefix.txt_sort_for_noSingles.bed > readsKeep.bed
rm readsPlusPrefix.txt_sort_for_noSingles.bed 
awk 'BEGIN{OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' readsKeep.bed > ${sample}_reads.onlyPairs.bed
rm readsKeep.bed

# make sure reasonable read numbers
echo "Number of reads in top1000 bed file, after removing singles:"
echo $(wc -l ${sample}_reads.onlyPairs.bed)

# split up file to run MM analysis
mkdir ${sample}_MM_temp_turnover
split --lines=10000 ${sample}_reads.onlyPairs.bed ${sample}_MM_temp_turnover/tmp.

# NOW RUN AGAIN FOR BOTTOM 500 GENES

timepoint=5

# run bottom500genes_turnover.R
/n/groups/churchman/bms36/programs/bottom500genes_turnover.R $gs_out $ptc_list $gtf_list $sample $timepoint

# get reads that align to the bottom500 genes
samtools view -@ 3 -b -L ${sample}_bottom500genes_turnover.bed ${sample}_sort.human_noMT.bam > ${sample}_bottom500genes_turnover.bam

#echo "Running samtools flagstat: all input reads"
samtools flagstat -@ 3 ${sample}_sort.human_noMT.bam

#echo "Running samtools flagstat: top 100 genes"
samtools flagstat -@ 3 ${sample}_bottom500genes_turnover.bam

# convert from .bam to .sam
samtools view -@ 3 -h -O SAM -o ${sample}_bottom500genes.paired.sam ${sample}_bottom500genes_turnover.bam
rm ${sample}_bottom500genes_turnover.bam

# remove soft clipped bases from reads
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/biostar84452.jar -o ${sample}_bottom500genes.paired.noSoft.sam ${sample}_bottom500genes.paired.sam
rm ${sample}_bottom500genes.paired.sam

# converting cigar string to discriminate between mismatches and matches
java -jar /n/groups/churchman/bms36/programs/jvarkit/dist/samfixcigar.jar -r $genome_fasta ${sample}_bottom500genes.paired.noSoft.sam > ${sample}_bottom500genes.paired.noSoft.fixCigar.sam
rm ${sample}_bottom500genes.paired.noSoft.sam

# get rid of sam header
samtools view -@ 3 -O SAM -o ${sample}_bottom500genes.paired.noSoft.fixCigar.noHead.sam ${sample}_bottom500genes.paired.noSoft.fixCigar.sam

# convert to bam (the sam with header)
samtools view -@ 3 -b -o ${sample}_bottom500genes.paired.noSoft.fixCigar.bam ${sample}_bottom500genes.paired.noSoft.fixCigar.sam
rm ${sample}_bottom500genes.paired.noSoft.fixCigar.sam

# convert to modified bed format
bedtools bamtobed -cigar -i ${sample}_bottom500genes.paired.noSoft.fixCigar.bam > ${sample}_reads_bottom500genes.bed
rm ${sample}_bottom500genes.paired.noSoft.fixCigar.bam

# run modifyBed.R script to edit bed file so MM analysis can be run on it *Need more than 25G memory for this step
/n/groups/churchman/bms36/programs/modifyBed.R ${sample}_reads_bottom500genes.bed ${sample}_bottom500genes.paired.noSoft.fixCigar.noHead.sam
rm ${sample}_reads_bottom500genes.bed
rm ${sample}_bottom500genes.paired.noSoft.fixCigar.noHead.sam

# make sure the resulting file from modifyBed.R contains same read numbers as above
echo "Number of reads in bottom500 bed file, before removing singles:"
echo $(wc -l ${sample}_reads_bottom500genes.bedwithC.bed)

# keep only reads that are still paired correctly
cut -f 4 ${sample}_reads_bottom500genes.bedwithC.bed | awk -F'/' 'BEGIN{OFS="\t";} {print $1}' > readPrefix.txt
paste ${sample}_reads_bottom500genes.bedwithC.bed readPrefix.txt > readsPlusPrefix.txt
rm readPrefix.txt
rm ${sample}_reads_bottom500genes.bedwithC.bed
/n/groups/churchman/bms36/sandbox/TimelapseSeq/MismatchScripts_batchSubmission/sortReads_byName_noSingles.R readsPlusPrefix.txt
rm readsPlusPrefix.txt
uniq -f 10 -D readsPlusPrefix.txt_sort_for_noSingles.bed > readsKeep.bed
rm readsPlusPrefix.txt_sort_for_noSingles.bed 
awk 'BEGIN{OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' readsKeep.bed > ${sample}_reads.onlyPairs.bed
rm readsKeep.bed

# make sure reasonable read numbers
echo "Number of reads in bottom500 bed file, after removing singles:"
echo $(wc -l ${sample}_reads.onlyPairs.bed)

# split up file to run MM analysis
mkdir ${sample}_MM_temp_turnover_slow
split --lines=10000 ${sample}_reads.onlyPairs.bed ${sample}_MM_temp_turnover_slow/tmp.

