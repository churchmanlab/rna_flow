#!/bin/bash
#Robert Ietswaart

date=20240309
baseDir=/n/groups/churchman/bms36/2022-11-04_RBP_KDs_ii/
resourceDir=/n/groups/churchman/ri23/bseq/GS${date}_KD_ii/
# scDir=/n/groups/churchman/ri23/code/
# notifEm=robert_ietswaart@hms.harvard.edu
# load required modules
module load gcc/6.2.0
module load java/jdk-1.8u112

mkdir -p ${resourceDir}/LogErr

samples='DIS3' #'EXOSC10 PABPN1 ZFC3H1 scramble' #'scramble' # 
compartments='nuc tot' #'nuc' #
reps='rep2' #'rep1 rep2' #
read_types='exons exons_bamlist_v2.0.5d_oml unspliced_junctions introns retained_introns retained_introns_with_exons protein_coding' #not_introns no_retained_introns

#Estimates from GS20240305_KD_ii_median_difference_w_scr.xlsx
#Assumptions: 
#1. median MAP fraction new RNA over all expressed genes in K562 scramble best matches K562 WT: determines pc_scr (per rep)
#2. median MAP fraction new RNA over all expressed genes difference between KD and scr = 0: determines ratio pc_scr / pc_KD
pc_scr_r1='0.0082'
pc_scr_r2='0.0068'
pc_DIS3_r1='0.0137'
pc_DIS3_r2='0.0068'
pc_EXOSC10_r1='0.0137'
pc_EXOSC10_r2='0.017'
pc_PABPN1_r1='0.0103'
pc_PABPN1_r2='0.0113'
pc_ZFC3H1_r1='0.0103'
pc_ZFC3H1_r2='0.0113'

for x in $samples
do
    for rep in $reps
    do 
        for comp in $compartments
        do
            for rt in $read_types
            do
                target_folder=${baseDir}GS_2023/${x}/${comp}/${rep}/${date}_${rt}/
                echo $target_folder
                mkdir -p $target_folder
                cd $target_folder
                rm *

                bamlist=${x}_${comp}_${rep}_noMT_${rt}.bamlist

                if [ ${comp} == 'tot' ]
                then
                    if [[ ${rt} == 'unspliced_junctions' || ${rt} == 'introns' || ${rt} == 'not_introns' ]]
                    then
                        echo "${baseDir}STAR_2023/${x}_0_${comp}_${rep}/${rt}/${x}_UL_${comp}_${rep}_sort.human_noMT_${rt}.bam" > $bamlist
                        echo "${baseDir}STAR_2023/${x}_60_${comp}_${rep}/${rt}/${x}_60_${comp}_${rep}_sort.human_noMT_${rt}.bam" >> $bamlist
                    else
                        echo "${baseDir}STAR_2023/${x}_0_${comp}_${rep}/${x}_UL_${comp}_${rep}_sort.human_noMT.bam" > $bamlist
                        echo "${baseDir}STAR_2023/${x}_60_${comp}_${rep}/${x}_60_${comp}_${rep}_sort.human_noMT.bam" >> $bamlist
                    fi

                    #Write rates file with pc and p_err values
                    pc_prefix=${x}_${comp}_${rep}_${date}_pc
                    echo -e "Rate\t${x}_UL\t${x}_60" > ${pc_prefix}.tsv  
                    echo -e "single_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                    if [[ ${x} == 'scramble' && ${rep} == 'rep1' ]]
                    then
                        echo -e "single_new\t${pc_scr_r1}\t${pc_scr_r1}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_scr_r1}\t${pc_scr_r1}" >> ${pc_prefix}.tsv  
                    elif [[ ${x} == 'scramble' && ${rep} == 'rep2' ]]
                    then
                        echo -e "single_new\t${pc_scr_r2}\t${pc_scr_r2}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_scr_r2}\t${pc_scr_r2}" >> ${pc_prefix}.tsv  
                    elif [[ ${x} == 'DIS3' && ${rep} == 'rep1' ]]
                    then
                        echo -e "single_new\t${pc_DIS3_r1}\t${pc_DIS3_r1}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_DIS3_r1}\t${pc_DIS3_r1}" >> ${pc_prefix}.tsv
                    elif [[ ${x} == 'DIS3' && ${rep} == 'rep2' ]]
                    then
                        echo -e "single_new\t${pc_DIS3_r2}\t${pc_DIS3_r2}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_DIS3_r2}\t${pc_DIS3_r2}" >> ${pc_prefix}.tsv
                    elif [[ ${x} == 'EXOXC10' && ${rep} == 'rep1' ]]
                    then
                        echo -e "single_new\t${pc_EXOSC10_r1}\t${pc_EXOSC10_r1}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_EXOSC10_r1}\t${pc_EXOSC10_r1}" >> ${pc_prefix}.tsv 
                    elif [[ ${x} == 'EXOXC10' && ${rep} == 'rep2' ]]
                    then
                        echo -e "single_new\t${pc_EXOSC10_r2}\t${pc_EXOSC10_r2}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_EXOSC10_r2}\t${pc_EXOSC10_r2}" >> ${pc_prefix}.tsv 
                    elif [[ ${x} == 'PABPN1' && ${rep} == 'rep1' ]]
                    then
                        echo -e "single_new\t${pc_PABPN1_r1}\t${pc_PABPN1_r1}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_PABPN1_r1}\t${pc_PABPN1_r1}" >> ${pc_prefix}.tsv
                    elif [[ ${x} == 'PABPN1' && ${rep} == 'rep2' ]]
                    then
                        echo -e "single_new\t${pc_PABPN1_r2}\t${pc_PABPN1_r2}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_PABPN1_r2}\t${pc_PABPN1_r2}" >> ${pc_prefix}.tsv
                    elif [[ ${x} == 'ZFC3H1' && ${rep} == 'rep1' ]]
                    then
                        echo -e "single_new\t${pc_ZFC3H1_r1}\t${pc_ZFC3H1_r1}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_ZFC3H1_r1}\t${pc_ZFC3H1_r1}" >> ${pc_prefix}.tsv 
                    elif [[ ${x} == 'ZFC3H1' && ${rep} == 'rep2' ]]
                    then
                        echo -e "single_new\t${pc_ZFC3H1_r2}\t${pc_ZFC3H1_r2}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_ZFC3H1_r2}\t${pc_ZFC3H1_r2}" >> ${pc_prefix}.tsv 
                    fi #sample

                else #NUC RNA 
                    if [[ ${rt} == 'unspliced_junctions' || ${rt} == 'introns' || ${rt} == 'not_introns' ]]
                    then
                        echo "${baseDir}STAR_2023/${x}_60_${comp}_${rep}/${rt}/${x}_60_${comp}_${rep}_sort.human_noMT_${rt}.bam" > $bamlist
                    else
                        echo "${baseDir}STAR_2023/${x}_60_${comp}_${rep}/${x}_60_${comp}_${rep}_sort.human_noMT.bam" > $bamlist
                    fi
                    
                    #Write rates file with pc and p_err values
                    pc_prefix=${x}_${comp}_${rep}_${date}_pc
                    echo -e "Rate\t${x}_60_${comp}_${rep}_sort.human_noMT" > ${pc_prefix}.tsv
                    echo -e "single_old\t0.001" >> ${pc_prefix}.tsv
                    if [[ ${x} == 'scramble' && ${rep} == 'rep1' ]]
                    then
                        echo -e "single_new\t${pc_scr_r1}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_scr_r1}" >> ${pc_prefix}.tsv  
                    elif [[ ${x} == 'scramble' && ${rep} == 'rep2' ]]
                    then
                        echo -e "single_new\t${pc_scr_r2}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_scr_r2}" >> ${pc_prefix}.tsv  
                    elif [[ ${x} == 'DIS3' && ${rep} == 'rep1' ]]
                    then
                        echo -e "single_new\t${pc_DIS3_r1}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_DIS3_r1}" >> ${pc_prefix}.tsv 
                    elif [[ ${x} == 'DIS3' && ${rep} == 'rep2' ]]
                    then
                        echo -e "single_new\t${pc_DIS3_r2}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_DIS3_r2}" >> ${pc_prefix}.tsv 
                    elif [[ ${x} == 'EXOSC10' && ${rep} == 'rep1' ]]
                    then
                        echo -e "single_new\t${pc_EXOSC10_r1}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_EXOSC10_r1}" >> ${pc_prefix}.tsv
                    elif [[ ${x} == 'EXOSC10' && ${rep} == 'rep2' ]]
                    then
                        echo -e "single_new\t${pc_EXOSC10_r2}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_EXOSC10_r2}" >> ${pc_prefix}.tsv
                    elif [[ ${x} == 'PABPN1' && ${rep} == 'rep1' ]]
                    then
                        echo -e "single_new\t${pc_PABPN1_r1}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_PABPN1_r1}" >> ${pc_prefix}.tsv  
                    elif [[ ${x} == 'PABPN1' && ${rep} == 'rep2' ]]
                    then
                        echo -e "single_new\t${pc_PABPN1_r2}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_PABPN1_r2}" >> ${pc_prefix}.tsv  
                    elif [[ ${x} == 'ZFC3H1' && ${rep} == 'rep1' ]]
                    then
                        echo -e "single_new\t${pc_ZFC3H1_r1}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_ZFC3H1_r1}" >> ${pc_prefix}.tsv  
                    elif [[ ${x} == 'ZFC3H1' && ${rep} == 'rep2' ]]
                    then
                        echo -e "single_new\t${pc_ZFC3H1_r2}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc_ZFC3H1_r2}" >> ${pc_prefix}.tsv      
                      
                    fi #sample
                fi #comp
            done
        done
    done
done
