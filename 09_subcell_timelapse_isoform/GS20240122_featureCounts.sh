#!/bin/bash
#Robert Ietswaart 
#generate a gene coverage file using featureCounts
#Source: NucWash_featureCounts_20231005.sh
#Run in interactive session

date=20240216
scDir=/n/groups/churchman/ri23/code/bseq/
cd ${scDir}

baseDir=/n/groups/churchman/ri23/bseq/GS20240122_KD_ii/
initDir=/n/groups/churchman/bms36/2022-11-04_RBP_KDs_ii/STAR_2023/
notifEm=robert_ietswaart@hms.harvard.edu

program=featureCounts

# load modules
module load gcc/6.2.0
module load subread/1.6.2

nthread=1 #slurm: 4

outDir=${baseDir}${program}/
mkdir -p ${outDir}LogErr

samples='DIS3 EXOSC10 PABPN1 ZFC3H1 scramble'
reps='rep1 rep2'
read_types='all unspliced_junctions introns not_introns'
out_types="gene exon alt"

bam_files1=""
bam_files2=""
for rt in $read_types
do
    for x in $samples
    do
        for rep in $reps
        do 
            if [ $rt == 'all' ]
            then
                # TOTAL RNA
                bam_files1="${bam_files1} ${initDir}${x}_0_tot_${rep}/${x}_UL_tot_${rep}_sort.human_noMT.bam"
                bam_files1="${bam_files1} ${initDir}${x}_60_tot_${rep}/${x}_60_tot_${rep}_sort.human_noMT.bam"

                #NUC RNA  
                bam_files1="${bam_files1} ${initDir}${x}_60_nuc_${rep}/${x}_60_nuc_${rep}_sort.human_noMT.bam"
            else
                # TOTAL RNA
                bam_files2="${bam_files2} ${initDir}${x}_0_tot_${rep}/${rt}/${x}_UL_tot_${rep}_sort.human_noMT_${rt}.bam"
                bam_files2="${bam_files2} ${initDir}${x}_60_tot_${rep}/${rt}/${x}_60_tot_${rep}_sort.human_noMT_${rt}.bam"

                #NUC RNA  
                bam_files2="${bam_files2} ${initDir}${x}_60_nuc_${rep}/${rt}/${x}_60_nuc_${rep}_sort.human_noMT_${rt}.bam"
            fi
        done
    done
done

echo $bam_files1

echo $bam_files2

#By default, featureCounts does not count reads overlapping with more than one feature 
#(or more than one meta-feature when summarizing at meta-feature level). 
#Users can use the -O option to instruct featureCounts to count such reads 
#(they will be assigned to all their overlapping features or meta-features). 
#--fracOverlap Minimum fraction of overlapping bases in a read that is required for read assignment. 
#Value should be within range [0,1]. 0 by default. Use 1 for exon to mimick GRAND-SLAM results

#for exon and gene level quantification, need regular gtf
gtfFile=/n/groups/churchman/ri23/bseq/GS20240126/K562_ensGRCh38_MTmod_dm6_ercc_cat.gtf
FEATURETYPE="gene exon"
for ft in $FEATURETYPE
do
    ${program} -a ${gtfFile} -g gene_id -t ${ft} -T ${nthread} -s 1 -d 20 --fracOverlap 1 -o ${outDir}GS20240122_KD_ii_${program}_${ft}.tsv ${bam_files1} 2> ${outDir}LogErr/featureCounts_${date}_${ft}.log;
done

#for the alternative (alt) read_types: need custom gtf
gtfFile=/n/groups/churchman/ri23/bseq/GS20240126/K562_ensGRCh38_MTmod_dm6_ercc_cat_genes_only.gtf
ft=exon #because custom gtf used labels the gene as "exon", needed for GRAND-SLAM
${program} -a ${gtfFile} -g gene_id -t ${ft} -T ${nthread} -s 1 -d 20 -o ${outDir}GS20240122_KD_ii_${program}_alt.tsv ${bam_files2} 2> ${outDir}LogErr/featureCounts_${date}_alt.log; 

#clean up featureCounts output files to count matrix format with headers corresponding to sample IDs.
for ot in $out_types 
do
    cut -f1,7- ${outDir}GS20240122_KD_ii_${program}_${ot}.tsv > ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;
    sed -i '1d' ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;
    sed -i '1s/.bam//g' ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;
    sed -i '1s/_sort.human_noMT//g' ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;
    initDir_with_slashes=$(echo ${initDir} | sed 's/\//\\\//g');
    sed -i "1s/${initDir_with_slashes}//g" ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;

    sed '1s/.bam//g' ${outDir}GS20240122_KD_ii_${program}_${ot}.tsv.summary > ${outDir}GS20240122_KD_ii_${program}_${ot}_summary.tsv;
    sed -i '1s/_sort.human_noMT//g' ${outDir}GS20240122_KD_ii_${program}_${ot}_summary.tsv;
    sed -i "1s/${initDir_with_slashes}//g" ${outDir}GS20240122_KD_ii_${program}_${ot}_summary.tsv;

    for rt in $read_types
    do
        for x in $samples
        do
            for rep in $reps
            do 
                if [ $ot == 'alt' ]
                then
                
                    #TOTAL RNA
                    sed -i "1s/${x}_0_tot_${rep}\/${rt}\///g" ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;
                    sed -i "1s/${x}_0_tot_${rep}\/${rt}\///g" ${outDir}GS20240122_KD_ii_${program}_${ot}_summary.tsv; 
                    sed -i "1s/${x}_60_tot_${rep}\/${rt}\///g" ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;
                    sed -i "1s/${x}_60_tot_${rep}\/${rt}\///g" ${outDir}GS20240122_KD_ii_${program}_${ot}_summary.tsv; 

                    #NUC RNA                 
                    sed -i "1s/${x}_60_nuc_${rep}\/${rt}\///g" ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;
                    sed -i "1s/${x}_60_nuc_${rep}\/${rt}\///g" ${outDir}GS20240122_KD_ii_${program}_${ot}_summary.tsv;  
                
                else               
                
                    #TOTAL RNA
                    sed -i "1s/${x}_0_tot_${rep}\/${x}_UL_tot_${rep}/${x}_UL_tot_${rep}_${ot}/g" ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;
                    sed -i "1s/${x}_0_tot_${rep}\/${x}_UL_tot_${rep}/${x}_UL_tot_${rep}_${ot}/g" ${outDir}GS20240122_KD_ii_${program}_${ot}_summary.tsv;
                    sed -i "1s/${x}_60_tot_${rep}\/${x}_60_tot_${rep}/${x}_60_tot_${rep}_${ot}/g" ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;
                    sed -i "1s/${x}_60_tot_${rep}\/${x}_60_tot_${rep}/${x}_60_tot_${rep}_${ot}/g" ${outDir}GS20240122_KD_ii_${program}_${ot}_summary.tsv;

                    #NUC
                    sed -i "1s/${x}_60_nuc_${rep}\/${x}_60_nuc_${rep}/${x}_60_nuc_${rep}_${ot}/g" ${outDir}GS20240122_KD_ii_counts_${ot}.tsv;
                    sed -i "1s/${x}_60_nuc_${rep}\/${x}_60_nuc_${rep}/${x}_60_nuc_${rep}_${ot}/g" ${outDir}GS20240122_KD_ii_${program}_${ot}_summary.tsv;              

                fi
            done
        done
    done
done

echo 'end of featureCounts sh script'
