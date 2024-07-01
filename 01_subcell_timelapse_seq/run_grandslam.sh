#!/bin/bash

#SBATCH -c 4
#SBATCH -t 720
#SBATCH -p short
#SBATCH --mem=25G
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=smalec@g.harvard.edu

# see: https://github.com/erhard-lab/gedi/wiki/GRAND-SLAM

# https://www.onwardpath.com/2018/07/09/1368.html

module load gcc/6.2.0
module load java/jdk-1.8u112
module load R/3.5.1

sample=`pwd | cut -d '/' -f 8`
prefix=`ls *.snpdata | cut -d '.' -f 1`

/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -prefix $prefix -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ../${sample}_noMT.cit
