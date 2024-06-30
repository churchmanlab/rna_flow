#!/bin/bash
#Robert Ietswaart 
#sbatch job to generate STAR index of hg38 with yeast spike ins 

date=20230915

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/
baseDir=/n/groups/churchman/ri23/genomes/

#1. Concatenate human / yeast genome fasta file
#fasta for 11 yeast spike-in transcripts received from Karine Choquet. 
#each transcript is its own contig/chromosome: sacCer3/sacCer3_pool1.fa
cd ${baseDir}
cat hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa sacCer3/sacCer3_pool1.fa >> hg38sacCer3/Homo_sapiens.GRCh38.dna.primary_assembly.and.sacCer3_pool1.fa

# 2. Concatenate human genome / yeast spike-in gtf files
# Note the yeast gtf is made manually corresponding to the yeast fasta file
cd ${baseDir}
filename_hg38=Homo_sapiens.GRCh38.86.gtf
filename_sacCer3=sacCer3_pool1.gtf
output_file=GRCh38.86_sacCer3_pool1.gtf
cat hg38/$filename_hg38 sacCer3/$filename_sacCer3 >> hg38sacCer3/$output_file

#3. Create STAR index  
mkdir -p ${baseDir}/hg38sacCer3/STARindex/LogErr
outDir=${baseDir}/hg38sacCer3/STARindex/
fastaFile=${baseDir}/hg38sacCer3/Homo_sapiens.GRCh38.dna.primary_assembly.and.sacCer3_pool1.fa
gtfFile=${baseDir}/hg38sacCer3/GRCh38.86_sacCer3_pool1.gtf
program=STAR
module load gcc/6.2.0
module load star/2.7.0a
nthread=4
sbatch -n ${nthread} --mem-per-cpu=16G -t 0-03:00:00 --job-name=hg38sacCer3_STARindex \
   -o ${outDir}/LogErr/hg38sacCer3_STARindex.log -e ${outDir}/LogErr/hg38sacCer3_STARindex.err  \
   -p short --mail-user=${notifEm} --mail-type=ALL \
   --wrap="${program} --runMode genomeGenerate --genomeDir ${outDir} \
           --genomeFastaFiles ${fastaFile} --sjdbGTFfile ${gtfFile} \
           --runThreadN ${nthread}"
