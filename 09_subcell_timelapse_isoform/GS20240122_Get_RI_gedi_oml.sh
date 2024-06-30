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

fa=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_dm6_ercc_cat.fasta

read_types='default protein_coding retained_introns retained_introns_with_exons no_retained_introns' #'retained_introns' #

for rt in $read_types
do
    #Source gtfs: GS20240122_Get_retained_intron_gtf.sh
    if [ ${rt} == 'default' ]
    then
        gtf=K562_ensGRCh38_MTmod_dm6_ercc_cat
    else
        gtf=K562_ensGRCh38_MTmod_dm6_ercc_cat_${rt}   
    fi

    #Generate a GRAND-SLAM index file (.oml) format
    echo "begin ${rt}"
    $program -e IndexGenome -s $fa -a ${gtf}.gtf -nobowtie -nokallisto -nostar -p -D 
    #unnecessary arguments
    #-D output debugging info
    #-p output progress

    #This generates the .oml file in home 
    cp /home/ri23/.gedi/genomic/${gtf}.oml ./

done

echo "end of script"