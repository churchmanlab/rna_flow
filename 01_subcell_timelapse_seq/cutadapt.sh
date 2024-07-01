#!/bin/bash

#SBATCH -c 4
#SBATCH -t 1440
#SBATCH -p priority
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=smalec@g.harvard.edu

### LAST UPDATED:		04/06/2020 by Brendan Smalec

# load modules
module load gcc/6.2.0
module load python/3.7.4
module load cutadapt/2.5

files='U5'

# set directory for raw files (no trailing slash)
raw_dir='/n/groups/churchman/bms36/2021-04-07_T_U/raw_reads/nuc_tot/210402_A00794_0397_BH2GMLDRXY_i5revcomp/fastq/2021-03-26_BMS'

adapter_R1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
adapter_R2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

for x in $files

do 

# tot

# Cutadapt to trim adapters, allowing 20% error rate since read quality is low at 3' end, require minimum length 20 (-m 20), trim 3 nt from 5' end, trim Gs that result from sequence runoff [no color, black] (--nextseq-trim=5), max 0 Ns allowed in sequence (max-n 0)
cutadapt -a ${adapter_R1} -A ${adapter_R2} -e .2 -m 20 -u 3 -j 4 -o ${x}_tot_R1_noadapt.fastq.gz -p ${x}_tot_R2_noadapt.fastq.gz ${raw_dir}/${x}-tot-R1.fastq.gz ${raw_dir}/${x}-tot-R2.fastq.gz

# Quality trim at both ends, but adding this quality filter to the above command causes quality trimming before adaptor trimming.
cutadapt -a ${adapter_R1} -A ${adapter_R2} -q 20,0 --nextseq-trim=20 -m 20 -j 4 -o ${x}_tot_R1_trimmed.fastq.gz -p ${x}_tot_R2_trimmed.fastq.gz ${x}_tot_R1_noadapt.fastq.gz ${x}_tot_R2_noadapt.fastq.gz

# remove intermediates
rm ${x}_tot_R1_noadapt.fastq.gz
rm ${x}_tot_R2_noadapt.fastq.gz

# And finally remove 3nt from 3' end of Read2 only to match the 3 trimmed from the 5' end of Read1
cutadapt -j 4 -m 20 -U -3 --max-n 0 -o ${x}_tot_R1_cutadapt.fastq.gz -p ${x}_tot_R2_cutadapt.fastq.gz ${x}_tot_R1_trimmed.fastq.gz ${x}_tot_R2_trimmed.fastq.gz

# remove intermediates
rm ${x}_tot_R1_trimmed.fastq.gz
rm ${x}_tot_R2_trimmed.fastq.gz

# nuc

# Cutadapt to trim adapters, allowing 20% error rate since read quality is low at 3' end, require minimum length 20 (-m 20), trim 3 nt from 5' end, trim Gs that result from sequence runoff [no color, black] (--nextseq-trim=5), max 0 Ns allowed in sequence (max-n 0)
cutadapt -a ${adapter_R1} -A ${adapter_R2} -e .2 -m 20 -u 3 -j 4 -o ${x}_nuc_R1_noadapt.fastq.gz -p ${x}_nuc_R2_noadapt.fastq.gz ${raw_dir}/${x}-nuc-R1.fastq.gz ${raw_dir}/${x}-nuc-R2.fastq.gz

# Quality trim at both ends, but adding this quality filter to the above command causes quality trimming before adaptor trimming.
cutadapt -a ${adapter_R1} -A ${adapter_R2} -q 20,0 --nextseq-trim=20 -m 20 -j 4 -o ${x}_nuc_R1_trimmed.fastq.gz -p ${x}_nuc_R2_trimmed.fastq.gz ${x}_nuc_R1_noadapt.fastq.gz ${x}_nuc_R2_noadapt.fastq.gz

# remove intermediates
rm ${x}_nuc_R1_noadapt.fastq.gz
rm ${x}_nuc_R2_noadapt.fastq.gz

# And finally remove 3nt from 3' end of Read2 only to match the 3 trimmed from the 5' end of Read1
cutadapt -j 4 -m 20 -U -3 --max-n 0 -o ${x}_nuc_R1_cutadapt.fastq.gz -p ${x}_nuc_R2_cutadapt.fastq.gz ${x}_nuc_R1_trimmed.fastq.gz ${x}_nuc_R2_trimmed.fastq.gz

# remove intermediates
rm ${x}_nuc_R1_trimmed.fastq.gz
rm ${x}_nuc_R2_trimmed.fastq.gz

done