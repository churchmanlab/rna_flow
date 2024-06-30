#!/bin/bash
#Robert Ietswaart
#Source: GS20240122_KDii_run_GS.sh and GS20240122_KDii_RI_run_GS.sh

date=20240305
baseDir=/n/groups/churchman/bms36/2022-11-04_RBP_KDs_ii/
resourceDir=/n/groups/churchman/ri23/bseq/GS20240126/
# scDir=/n/groups/churchman/ri23/code/
notifEm=robert_ietswaart@hms.harvard.edu

nthread=4

# load required modules
module load gcc/6.2.0
module load java/jdk-1.8u112

samples='scramble' #'DIS3 EXOSC10 PABPN1 ZFC3H1 scramble' #'DIS3' #'scramble' #
compartments='nuc' #'nuc tot' #
reps='rep1' #'rep1 rep2' #
rt='exons_bamlist_v2.0.5d_oml'

PC_KDs='0.015 0.025 0.035'
RATIO_PC_SCR_KD='0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6'


for x in $samples
do
    for rep in $reps
    do
        for comp in $compartments
        do
            oml=${resourceDir}K562_ensGRCh38_MTmod_dm6_ercc_cat.oml
            reads_input=../${x}_${comp}_${rep}_noMT_${rt}.bamlist
            
            if [ ${x} == 'scramble' ]
            then
                for pc_kd in $PC_KDs
                do 
                    for r in $RATIO_PC_SCR_KD 
                    do
                        pc=$(echo "scale=10; $r * $pc_kd" | bc | sed 's/0*$//;s/\.$//')
                        [[ "$pc" == .* ]] && pc="0$pc" # Ensure leading zero if necessary    
                        pc_prefix=pc_${pc}_version_pc_kd_${pc_kd}
                        pc_input=../${pc_prefix}.tsv
                        prefix=${x}_${comp}_${rep}_${rt}_${pc_prefix}
                        
                        target_folder=${baseDir}GS_2023/${x}/${comp}/${rep}/${date}_${rt}/${pc_prefix}/    
                        cd $target_folder
                        rm *
                        
                        # submit job
                        OUT=run_${prefix}_GS.sh #THIS is the name of the script to be written.
                        echo "#!/bin/bash" > $OUT
                        echo "#SBATCH -c ${nthread}" >> $OUT
                        echo "#SBATCH -t 300" >> $OUT #300
                        echo "#SBATCH -p short" >> $OUT
                        echo "#SBATCH --mem=16G" >> $OUT #8 leads to OUT_OF_MEMORY for scramble nuc rep 1
                        echo "#SBATCH -o %j.out" >> $OUT
                        echo "#SBATCH -e %j.err" >> $OUT
                        echo "#SBATCH --mail-type=ALL" >> $OUT
                        echo "#SBATCH --mail-user=${notifEm}" >> $OUT
                        echo "" >> $OUT
                        echo "module load gcc/6.2.0" >> $OUT
                        echo "module load java/jdk-1.8u112" >> $OUT
                        echo "module load R/3.5.1" >> $OUT
                        echo "" >> $OUT
                        if [ $comp == 'nuc' ]
                        then
                            echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic ${oml} -prefix ${prefix} -mode All -snpConv 0.2 -nthreads ${nthread} -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ${reads_input}" >> $OUT
                        elif [ $comp == 'tot' ]
                        then
                            echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic ${oml} -prefix ${prefix} -no4sUpattern UL -mode All -snpConv 0.2 -nthreads ${nthread} -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ${reads_input}" >> $OUT
                        fi
                        echo "" >> $OUT
                        echo "cp ${pc_input} ${prefix}.rates.tsv" >> $OUT
                        echo "mv ${prefix}.tsv ${prefix}_run_i.tsv" >> $OUT
                        echo "rm ${prefix}.ext.tsv " >> $OUT
                        echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic ${oml} -prefix ${prefix} -mode All -snpConv 0.2 -nthreads ${nthread} -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ${reads_input}" >> $OUT      
                        echo "cp ${prefix}.tsv /n/groups/churchman/ri23/bseq/GS${date}_KD_ii/GS${date}_${prefix}.tsv" >> $OUT
                        chmod u+x $OUT  # MUST change permission, otherwise will not submit
                        sbatch $OUT  # Submit the script   
                        
                    done
                done   
            else #KD
                for pc in $PC_KDs
                do 
                
                    pc_prefix=pc_${pc}
                    pc_input=../${pc_prefix}.tsv
                    prefix=${x}_${comp}_${rep}_${rt}_${pc_prefix}
                    
                    target_folder=${baseDir}GS_2023/${x}/${comp}/${rep}/${date}_${rt}/${pc_prefix}/
                    cd $target_folder
                    rm *

                    # submit job
                    OUT=run_${prefix}_GS.sh #THIS is the name of the script to be written.
                    echo "#!/bin/bash" > $OUT
                    echo "#SBATCH -c ${nthread}" >> $OUT
                    echo "#SBATCH -t 300" >> $OUT #300
                    echo "#SBATCH -p short" >> $OUT
                    echo "#SBATCH --mem=16G" >> $OUT
                    echo "#SBATCH -o %j.out" >> $OUT
                    echo "#SBATCH -e %j.err" >> $OUT
                    echo "#SBATCH --mail-type=ALL" >> $OUT
                    echo "#SBATCH --mail-user=${notifEm}" >> $OUT
                    echo "" >> $OUT
                    echo "module load gcc/6.2.0" >> $OUT
                    echo "module load java/jdk-1.8u112" >> $OUT
                    echo "module load R/3.5.1" >> $OUT
                    echo "" >> $OUT
                    if [ $comp == 'nuc' ]
                    then
                        echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic ${oml} -prefix ${prefix} -mode All -snpConv 0.2 -nthreads ${nthread} -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ${reads_input}" >> $OUT
                    elif [ $comp == 'tot' ]
                    then
                        echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic ${oml} -prefix ${prefix} -no4sUpattern UL -mode All -snpConv 0.2 -nthreads ${nthread} -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ${reads_input}" >> $OUT
                    fi
                    echo "" >> $OUT
                    echo "cp ${pc_input} ${prefix}.rates.tsv" >> $OUT
                    echo "mv ${prefix}.tsv ${prefix}_run_i.tsv" >> $OUT
                    echo "rm ${prefix}.ext.tsv " >> $OUT
                    echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic ${oml} -prefix ${prefix} -mode All -snpConv 0.2 -nthreads ${nthread} -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ${reads_input}" >> $OUT      
                    echo "cp ${prefix}.tsv /n/groups/churchman/ri23/bseq/GS${date}_KD_ii/GS${date}_${prefix}.tsv" >> $OUT
                    chmod u+x $OUT  # MUST change permission, otherwise will not submit
                    sbatch $OUT  # Submit the script
                    
                done
            fi                           
        done
    done
done