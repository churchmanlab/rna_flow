#!/bin/bash

#SBATCH -c 1
#SBATCH -t 80
#SBATCH -p priority
#SBATCH --mem=100G
#SBATCH -o %j.out
#SBATCH -e %j.err



# load required modules
module load gcc/6.2.0
module load R/3.5.1
module load bedtools/2.27.1

# cat
cat rpf*fragments.bed > reads_MM.fragments.bed
#rm tmp*
rm slurm*

# sort
bedtools sort -i reads_MM.fragments.bed > reads_MM.fragments.sort.bed
rm reads_MM.fragments.bed

# reformat
awk -F'\t' 'BEGIN{OFS="\t"; print "chr\tstart\tend\tstrand\tseq\tfrag_length\ttot_mismatches\tTC_mismatches";} {print $1"\t"$2"\t"$3"\t"$6"\t"$9"\t"$21"\t"$18"\t"$19}' reads_MM.fragments.sort.bed > frag_MMfrequency_All.txt

# rename based on path (add sample name to front)
x=`pwd | cut -d '/' -f 8`
mv frag_MMfrequency_All.txt ./${x}_frag_MMfrequency_All.txt