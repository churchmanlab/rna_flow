#!/bin/bash

#Author: Robert Ietswaart 
date=20240126
baseDir=/n/groups/churchman/bms36/2022-11-04_RBP_KDs_ii/

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/bseq/
cd ${scDir}

read_types='unspliced_junctions introns not_introns' 

#1. Generate new folders for new bam files
prefixes='DIS3_0_tot_rep1 DIS3_0_tot_rep2 DIS3_60_nuc_rep1 DIS3_60_nuc_rep2 DIS3_60_tot_rep1 DIS3_60_tot_rep2 EXOSC10_0_tot_rep1 EXOSC10_0_tot_rep2 EXOSC10_60_nuc_rep1 EXOSC10_60_nuc_rep2 EXOSC10_60_tot_rep1 EXOSC10_60_tot_rep2 PABPN1_0_tot_rep1 PABPN1_0_tot_rep2 PABPN1_60_nuc_rep1 PABPN1_60_nuc_rep2 PABPN1_60_tot_rep1 PABPN1_60_tot_rep2 ZFC3H1_0_tot_rep1 ZFC3H1_0_tot_rep2 ZFC3H1_60_nuc_rep1 ZFC3H1_60_nuc_rep2 ZFC3H1_60_tot_rep1 ZFC3H1_60_tot_rep2 scramble_0_tot_rep1 scramble_0_tot_rep2 scramble_60_nuc_rep1 scramble_60_nuc_rep2 scramble_60_tot_rep1 scramble_60_tot_rep2'

for x in $prefixes 
do 
    echo $x
    for rt in $read_types
    do
        mkdir ${baseDir}STAR_2023/${x}/${rt}
    done
done


#2. Generate new folders for GS output
samples='DIS3 EXOSC10 PABPN1 ZFC3H1 scramble'
compartments='nuc tot'
reps='rep1 rep2'
tb_types='top bottom'
for x in $samples
do
    for comp in $compartments
    do
        for rep in $reps
        do
            for rt in $read_types
            do
                for tb in $tb_types
                do
                    target_folder=${baseDir}/GS_2023/${x}/${comp}/${rep}/${rt}/${tb}/
                    echo $target_folder
                    mkdir -p $target_folder
#                     #OLD, DOES NOT WORK: Copy previous GS over 
#                     cp ${baseDir}/GS_2023/${x}/${comp}/${rep}/${tb}/${x}_${comp}_${rep}.* $target_folder
#                     #Remove the previous GS output file that was copied over
#                     rm ${target_folder}${x}_${comp}_${rep}.tsv
                done
            done
        done
    done
done
