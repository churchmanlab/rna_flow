#!/bin/bash
#Robert Ietswaart
#Source: /n/groups/churchman/bms36/2022-11-04_RBP_KDs_ii/GS_2023/run_submit_gs_final.sh

date=20240126
baseDir=/n/groups/churchman/bms36/2022-11-04_RBP_KDs_ii/
resourceDir=/n/groups/churchman/ri23/bseq/GS${date}/
# scDir=/n/groups/churchman/ri23/code/
notifEm=robert_ietswaart@hms.harvard.edu

nthread=4

# load required modules
module load gcc/6.2.0
module load java/jdk-1.8u112

samples='DIS3 EXOSC10 PABPN1 ZFC3H1 scramble' #'DIS3' #
compartments='nuc tot' #'nuc' #
reps='rep1 rep2' #'rep1' #
tb_types='top bottom' #'top' #

read_types='unspliced_junctions introns not_introns'

for x in $samples
do
    for rep in $reps
    do
        for comp in $compartments
        do
            for rt in $read_types
            do
                for tb in $tb_types
                do
                    target_folder=${baseDir}GS_2023/${x}/${comp}/${rep}/${rt}/${tb}/
                    cd $target_folder
                    rm *
                    bamlist=../${x}_${comp}_${rep}_noMT_${rt}.bamlist

                    old_prefix=${x}_${comp}_${rep}
                    prefix=${x}_${comp}_${rep}_${rt}_${tb}                  
                 
                    # submit job
                    OUT=run_${x}_${comp}_${rep}_${rt}_${tb}_GS.sh #THIS is the name of the script to be written.
                    echo "#!/bin/bash" > $OUT
                    echo "#SBATCH -c ${nthread}" >> $OUT
                    echo "#SBATCH -t 300" >> $OUT
                    echo "#SBATCH -p short" >> $OUT
                    echo "#SBATCH --mem=25G" >> $OUT
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
                        echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic ${resourceDir}K562_ensGRCh38_MTmod_dm6_ercc_cat_genes_only.oml -prefix ${prefix} -mode All -snpConv 0.2 -nthreads ${nthread} -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ${bamlist}" >> $OUT
                    elif [ $comp == 'tot'  ]
                    then
                        echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic ${resourceDir}K562_ensGRCh38_MTmod_dm6_ercc_cat_genes_only.oml -prefix ${prefix} -no4sUpattern UL -mode All -snpConv 0.2 -nthreads ${nthread} -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ${bamlist}" >> $OUT
                    fi
                    echo "" >> $OUT
                    echo "cp ${baseDir}GS_2023/${x}/${comp}/${rep}/${tb}/${old_prefix}.rates.tsv ${prefix}.rates.tsv" >> $OUT
                    echo "mv ${prefix}.tsv ${prefix}_run_i.tsv" >> $OUT #TEMP FOR TEST
#                     echo "rm ${prefix}.tsv" >> $OUT
                    echo "rm ${prefix}.ext.tsv " >> $OUT
                    echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic ${resourceDir}K562_ensGRCh38_MTmod_dm6_ercc_cat_genes_only.oml -prefix ${prefix} -mode All -snpConv 0.2 -nthreads ${nthread} -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ${bamlist}" >> $OUT                    
                    chmod u+x $OUT  # MUST change permission, otherwise will not submit
                    sbatch $OUT  # Submit the script 
                                       
                done
            done
        done
    done
done