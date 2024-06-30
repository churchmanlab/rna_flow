#!/bin/bash
#Robert Ietswaart 
#generate a gene coverage file using featureCounts
#Run in interactive session after NucWash_STAR_alignment_20231005.sh

date=20231012v3
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

program=featureCounts

# load modules
module load gcc/6.2.0
module load subread/1.6.2

nthread=1 #slurm: 4

gtfFile=/n/groups/churchman/ri23/genomes/hg38sacCer3/GRCh38.86_sacCer3_pool1.gtf

FEATURETYPE="gene exon"

bam_files=""
FILES=`echo -e $samples`
for f in $FILES
do
    bam_files="${bam_files} ${baseDir}STAR/${f}/${f}_sorted_indexed.bam"
done

echo $bam_files

#By default, featureCounts does not count reads overlapping with more than one feature 
#(or more than one meta-feature when summarizing at meta-feature level). 
#Users can use the -O option to instruct featureCounts to count such reads 
#(they will be assigned to all their overlapping features or meta-features). 
    
#clean up featureCounts output files to count matrix format with headers corresponding to sample IDs.


for ft in $FEATURETYPE
do
    echo $ft
    outDir=${baseDir}${program}_${ft}/
    mkdir -p ${outDir}LogErr

    ${program} -a ${gtfFile} -g gene_id -t ${ft} -T ${nthread} -s 1 -d 20 -o ${outDir}CG20230908_${program}_${ft}.tsv ${bam_files} 2> ${outDir}LogErr/featureCounts_${date}_${ft}.log; #20231012v2/3

    cut -f1,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 ${outDir}CG20230908_${program}_${ft}.tsv > ${outDir}CG20230908_counts_${ft}.tsv;
    sed -i '1d' ${outDir}CG20230908_counts_${ft}.tsv;
    sed -i '1s/_sorted_indexed.bam//g' ${outDir}CG20230908_counts_${ft}.tsv;
    baseDir_with_slashes=$(echo ${baseDir} | sed 's/\//\\\//g');
    sed -i "1s/${baseDir_with_slashes}STAR\///g" ${outDir}CG20230908_counts_${ft}.tsv;
    for f in $FILES; do sed -i "1s/${f}\///g" ${outDir}CG20230908_counts_${ft}.tsv; done

    sed '1s/_sorted_indexed.bam//g' ${outDir}CG20230908_${program}_${ft}.tsv.summary > ${outDir}CG20230908_${program}_${ft}_summary.tsv;
    baseDir_with_slashes=$(echo ${baseDir} | sed 's/\//\\\//g');
    sed -i "1s/${baseDir_with_slashes}STAR\///g" ${outDir}CG20230908_${program}_${ft}_summary.tsv;
    for f in $FILES; do sed -i "1s/${f}\///g" ${outDir}CG20230908_${program}_${ft}_summary.tsv; done

done

echo 'end of featureCounts sh script'


    #old:
    # # works direct in terminal,so without slurm:
    #${program} -a ${gtfFile} -g gene_id -T ${nthread} -s 1 -o ${outDir}CG20230908_${program}.tsv ${bam_files}; #20231006v1
    # ${program} -a ${gtfFile} -g gene_id -T ${nthread} --primary -O -s 1 -d 20 -o ${outDir}CG20230908_${program}.tsv ${bam_files} 2> ${outDir}LogErr/featureCounts_${date}.log; #20231009v1 / 
    # ${program} -a ${gtfFile} -g gene_id -t gene -T ${nthread} --primary -O -s 1 -d 20 -o ${outDir}CG20230908_${program}.tsv ${bam_files} 2> ${outDir}LogErr/featureCounts_${date}.log; #20231012v1: control to see if more nuclear reads are mapped with gene (exon+intron) than exon only  (20231009v1)

#last line does not work properly, giving up on this for now because above in terminal works 
#sbatch -n ${nthread} --mem=16G -t 0-00:15:00 --job-name=${program}${date} \
# -o ${outDir}LogErr/${program}${date}.log -e ${outDir}LogErr/${program}${date}.err \
# -p short --mail-user=${notifEm} --mail-type=ALL \
# --wrap="${program} -a ${gtfFile} -g gene_id -T ${nthread} -s 1 -o ${outDir}CG20230908_${program}.tsv ${bam_files}; \
# cut -f1,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25 ${outDir}CG20230908_${program}.tsv > ${outDir}CG20230908_counts.tsv; \
# sed -i '1d' ${outDir}CG20230908_counts.tsv; \
# sed -i '1s/_sorted_indexed.bam//g' ${outDir}CG20230908_counts.tsv; \
# baseDir_with_slashes=\$(echo ${baseDir} | sed 's/\//\\\\\//g'); \
# sed -i \"1s/\${baseDir_with_slashes}STAR\///g\" ${outDir}CG20230908_counts.tsv; \
# for f in $FILES; do sed -i \"1s\/\$f\\/\\//g\" ${outDir}CG20230908_counts.tsv; done"



