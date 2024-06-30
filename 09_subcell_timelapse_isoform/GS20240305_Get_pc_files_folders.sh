#!/bin/bash
#Robert Ietswaart

date=20240305
baseDir=/n/groups/churchman/bms36/2022-11-04_RBP_KDs_ii/
resourceDir=/n/groups/churchman/ri23/bseq/GS${date}_KD_ii/
# scDir=/n/groups/churchman/ri23/code/
# notifEm=robert_ietswaart@hms.harvard.edu
# load required modules
module load gcc/6.2.0
module load java/jdk-1.8u112

mkdir -p ${resourceDir}/LogErr

samples='DIS3 EXOSC10 PABPN1 ZFC3H1 scramble' #'scramble' # 
compartments='nuc tot' #'nuc' #
reps='rep1 rep2' #'rep1' # 
rt='exons_bamlist_v2.0.5d_oml'

PC_KDs='0.015 0.025 0.035'
RATIO_PC_SCR_KD='0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6'

for x in $samples
do
    for rep in $reps
    do 
        for comp in $compartments
        do
            target_folder=${baseDir}GS_2023/${x}/${comp}/${rep}/${date}_${rt}/
            echo $target_folder
            mkdir -p $target_folder
            cd $target_folder
        
            bamlist=${x}_${comp}_${rep}_noMT_${rt}.bamlist
            
            if [ ${comp} == 'tot' ]
            then
                echo "${baseDir}STAR_2023/${x}_0_${comp}_${rep}/${x}_UL_${comp}_${rep}_sort.human_noMT.bam" > $bamlist
                echo "${baseDir}STAR_2023/${x}_60_${comp}_${rep}/${x}_60_${comp}_${rep}_sort.human_noMT.bam" >> $bamlist

                #Write rates file with pc and p_err values
                if [ ${x} == 'scramble' ]
                then
                    for pc_kd in $PC_KDs
                    do 
                        for r in $RATIO_PC_SCR_KD 
                        do
                            pc=$(echo "scale=10; $r * $pc_kd" | bc | sed 's/0*$//;s/\.$//')
                            [[ "$pc" == .* ]] && pc="0$pc" # Ensure leading zero if necessary
        
                            
                            pc_prefix=pc_${pc}_version_pc_kd_${pc_kd}
                            mkdir ${pc_prefix}
                            
                            echo -e "Rate\t${x}_UL\t${x}_60" > ${pc_prefix}.tsv
                            echo -e "single_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                            echo -e "single_new\t${pc}\t${pc}" >> ${pc_prefix}.tsv
                            echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                            echo -e "double_new\t${pc}\t${pc}" >> ${pc_prefix}.tsv  
                            
                        done
                    done
                        
                else #KD
                    for pc in $PC_KDs
                    do 
                        pc_prefix=pc_${pc}
                        mkdir ${pc_prefix}
                            
                        echo -e "Rate\t${x}_UL\t${x}_60" > ${pc_prefix}.tsv
                        echo -e "single_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "single_new\t${pc}\t${pc}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc}\t${pc}" >> ${pc_prefix}.tsv        
                    done
                fi #sample

            else #NUC RNA 
                echo "${baseDir}STAR_2023/${x}_60_${comp}_${rep}/${x}_60_${comp}_${rep}_sort.human_noMT.bam" > $bamlist

                #Write rates file with pc and p_err values
                if [ ${x} == 'scramble' ]
                then
                    for pc_kd in $PC_KDs
                    do 
                        for r in $RATIO_PC_SCR_KD 
                        do
                            pc=$(echo "scale=10; $r * $pc_kd" | bc | sed 's/0*$//;s/\.$//')
                            [[ "$pc" == .* ]] && pc="0$pc" # Ensure leading zero if necessary
                            
                            pc_prefix=pc_${pc}_version_pc_kd_${pc_kd}
                            mkdir ${pc_prefix}
                            
                            echo -e "Rate\t${x}_60_${comp}_${rep}_sort.human_noMT" > ${pc_prefix}.tsv
                            echo -e "single_old\t0.001" >> ${pc_prefix}.tsv
                            echo -e "single_new\t${pc}" >> ${pc_prefix}.tsv
                            echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                            echo -e "double_new\t${pc}" >> ${pc_prefix}.tsv  
                            
                        done
                    done
                        
                else #KD
                    for pc in $PC_KDs
                    do 
                        pc_prefix=pc_${pc}
                        mkdir ${pc_prefix}
                            
                        echo -e "Rate\t${x}_60_${comp}_${rep}_sort.human_noMT" > ${pc_prefix}.tsv
                        echo -e "single_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "single_new\t${pc}" >> ${pc_prefix}.tsv
                        echo -e "double_old\t0.001" >> ${pc_prefix}.tsv
                        echo -e "double_new\t${pc}" >> ${pc_prefix}.tsv       
                    done
                fi #sample
            fi #comp
        done
    done
done
