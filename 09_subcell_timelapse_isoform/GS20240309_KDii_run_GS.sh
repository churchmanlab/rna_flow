#!/bin/bash
#Robert Ietswaart
#Source: GS20240122_KDii_run_GS.sh and GS20240122_KDii_RI_run_GS.sh

date=20240309
baseDir=/n/groups/churchman/bms36/2022-11-04_RBP_KDs_ii/
resourceDir=/n/groups/churchman/ri23/bseq/GS20240126/
# scDir=/n/groups/churchman/ri23/code/
notifEm=robert_ietswaart@hms.harvard.edu

nthread=4

# load required modules
module load gcc/6.2.0
module load java/jdk-1.8u112

samples='DIS3' #'EXOSC10 PABPN1 ZFC3H1 scramble' # 'scramble' #
compartments='nuc tot' #'nuc'
reps='rep2' #'rep1 rep2' # 
read_types='exons exons_bamlist_v2.0.5d_oml unspliced_junctions introns retained_introns retained_introns_with_exons protein_coding' #not_introns no_retained_introns


for x in $samples
do
    for rep in $reps
    do
        for comp in $compartments
        do
            for rt in $read_types
            do
                if [ ${rt} == 'exons' ]
                then
                    oml=/home/mc348/.gedi/genomic/K562_dm6.oml
                    reads_input=../${x}_${comp}_${rep}_noMT.cit
                elif [ ${rt} == 'exons_bamlist_v2.0.5d_oml' ]
                then
                    oml=${resourceDir}K562_ensGRCh38_MTmod_dm6_ercc_cat.oml
                    reads_input=${x}_${comp}_${rep}_noMT_${rt}.bamlist
                elif [[ ${rt} == 'unspliced_junctions' || ${rt} == 'introns' || ${rt} == 'not_introns' ]]
                then
                    oml=${resourceDir}K562_ensGRCh38_MTmod_dm6_ercc_cat_genes_only.oml 
                    reads_input=${x}_${comp}_${rep}_noMT_${rt}.bamlist
                else
                    oml=${resourceDir}K562_ensGRCh38_MTmod_dm6_ercc_cat_${rt}.oml
                    reads_input=${x}_${comp}_${rep}_noMT_${rt}.bamlist
                fi
                
                target_folder=${baseDir}GS_2023/${x}/${comp}/${rep}/${date}_${rt}/
                cd $target_folder

                pc_prefix=${x}_${comp}_${rep}_${date}_pc
                prefix=${x}_${comp}_${rep}_${rt}                  

                # submit job
                OUT=run_${x}_${comp}_${rep}_${rt}_GS.sh #THIS is the name of the script to be written.
                echo "#!/bin/bash" > $OUT
                echo "#SBATCH -c ${nthread}" >> $OUT
                echo "#SBATCH -t 300" >> $OUT
                echo "#SBATCH -p short" >> $OUT
                echo "#SBATCH --mem=12G" >> $OUT
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
                echo "cp ${pc_prefix}.tsv ${prefix}.rates.tsv" >> $OUT
                echo "mv ${prefix}.tsv ${prefix}_run_i.tsv" >> $OUT
                echo "rm ${prefix}.ext.tsv " >> $OUT
                echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic ${oml} -prefix ${prefix} -mode All -snpConv 0.2 -nthreads ${nthread} -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads ${reads_input}" >> $OUT  
                echo "cp ${prefix}.tsv /n/groups/churchman/ri23/bseq/GS${date}_KD_ii/GS${date}_${prefix}.tsv" >> $OUT
                chmod u+x $OUT  # MUST change permission, otherwise will not submit
                sbatch $OUT  # Submit the script 
                                       
            done
        done
    done
done