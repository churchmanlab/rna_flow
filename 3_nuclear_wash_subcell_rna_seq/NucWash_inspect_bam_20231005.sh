#!/bin/bash
#Robert Ietswaart 
#sbatch job to generate a gene coverage file using featureCounts
#Run after NucWash_STAR_alignment_20231002.sh

date=20231005
scDir=/n/groups/churchman/ri23/code/bseq/
cd ${scDir}
suffix="" #default: "" #alternative (not used): _barcode-mismatches0
parameterFiles=NucWash_CG20230908${suffix}_parameters.in
echo $parameterFiles 

baseDir=`grep "Base Directory"            $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
samples=`grep "Sample Names"           $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
indexDir=`grep "Index Directory"            $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
initDir=`grep "Initial Files Directory"     $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
notifEm=`grep "Notification Email"       $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`

program=samtools

# load modules
module load gcc/6.2.0
module load ${program}/1.15.1

nthread=4

YEAST_TRANSCRIPTS="NM_001179372.1__TIM44 NM_001180740.3__NPL3 NM_001181339.1__YGR210C NM_001181986.1__ICT1 NM_001181013.1__ARO2 NM_001181842.1__HIF1 NM_001181160.1__IMO32 NM_001179305.1__ENO2 NM_001179373.1__YKE4 NM_001181805.1__HMS2 NM_001179170.1__BCD1"

FILES=`echo -e $samples`
FILES=1A_WC_S1
for f in $FILES
do
    echo ${f}
    outDir=${baseDir}STAR_old2/${f}
    mkdir -p ${outDir}/LogErr
    cd ${outDir}
    #Sort and index reads
    #Count number of 
    
#     sbatch -n ${nthread} --mem=16G -t 0-01:00:00 --job-name=${program}_${f} \
#     -o ${outDir}/LogErr/${program}_${f}.log -e ${outDir}/LogErr/${program}_${f}.err \
#     -p short --mail-user=${notifEm} --mail-type=ALL \
#     --wrap="

#NO LONGER NEEDED: added to STAR alignment script.
#     $program sort -@ ${nthread} -o ${f}_sorted_indexed.bam ${f}_Aligned.sortedByCoord.out.bam ; \
#     $program index -@ ${nthread} ${f}_sorted_indexed.bam ; \

#     $program view -@ ${nthread} ${f}_sorted_indexed.bam NM_001179372.1__TIM44 NM_001180740.3__NPL3 NM_001181339.1__YGR210C NM_001181986.1__ICT1 NM_001181013.1__ARO2 NM_001181842.1__HIF1 NM_001181160.1__IMO32 NM_001179305.1__ENO2 NM_001179373.1__YKE4 NM_001181805.1__HMS2 NM_001179170.1__BCD1 | wc -l #head # 

    $program view -@ ${nthread} ${f}_sorted_indexed.bam $YEAST_TRANSCRIPTS | wc -l #head # 
    
# 
    for yt in $YEAST_TRANSCRIPTS
    do
    echo $yt 
    $program view -@ ${nthread} ${f}_sorted_indexed.bam $yt | head # wc -l
    done
#   "  
done



echo 'end of inspect bam sh script'