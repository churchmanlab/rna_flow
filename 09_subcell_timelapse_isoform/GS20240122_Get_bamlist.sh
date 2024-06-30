#!/bin/bash
#Robert Ietswaart

date=20240126
baseDir=/n/groups/churchman/bms36/2022-11-04_RBP_KDs_ii/
resourceDir=/n/groups/churchman/ri23/bseq/GS${date}/
# scDir=/n/groups/churchman/ri23/code/
# notifEm=robert_ietswaart@hms.harvard.edu

# load required modules
module load gcc/6.2.0
module load java/jdk-1.8u112


samples='DIS3 EXOSC10 PABPN1 ZFC3H1 scramble'
reps='rep1 rep2'
read_types='unspliced_junctions introns not_introns' #'introns not_introns' #unspliced_junctions #'test'


for x in $samples
do
    for rep in $reps
    do 
        for rt in $read_types
        do
            # TOTAL RNA
            cd ${baseDir}GS_2023/${x}/tot/${rep}/${rt}/
  
            bamlist=${x}_tot_${rep}_noMT_${rt}.bamlist 
            
            echo "${baseDir}STAR_2023/${x}_0_tot_${rep}/${rt}/${x}_UL_tot_${rep}_sort.human_noMT_${rt}.bam" > $bamlist
            echo "${baseDir}STAR_2023/${x}_60_tot_${rep}/${rt}/${x}_60_tot_${rep}_sort.human_noMT_${rt}.bam" >> $bamlist


            #NUC RNA
            cd ${baseDir}GS_2023/${x}/nuc/${rep}/${rt}/
            
            bamlist=${x}_nuc_${rep}_noMT_${rt}.bamlist 
            
            echo "${baseDir}STAR_2023/${x}_60_nuc_${rep}/${rt}/${x}_60_nuc_${rep}_sort.human_noMT_${rt}.bam" > $bamlist 

        done
    done
done