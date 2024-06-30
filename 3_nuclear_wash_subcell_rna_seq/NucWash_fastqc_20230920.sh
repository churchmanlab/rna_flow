#!/bin/bash
#Robert Ietswaart 
#sbatch job to do fastqc on Nuc Wash fastq files
#Run after bcl2fastq_20230912.sh

date=20230920
scDir=/n/groups/churchman/ri23/code/bseq/
cd ${scDir}
suffix="" #default: "" #also tested: _barcode-mismatches0
parameterFiles=NucWash_CG20230908${suffix}_parameters.in
echo $parameterFiles 

baseDir=`grep "Base Directory"            $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
samples=`grep "Sample Names"           $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
indexDir=`grep "Index Directory"            $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
initDir=`grep "Initial Files Directory"     $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
notifEm=`grep "Notification Email"       $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`


program1=fastqc

module load gcc/6.2.0
module load ${program1}/0.11.9


nthread=1

mkdir -p ${baseDir}${program1}${suffix}/
FILES=`echo -e $samples`
for f in $FILES
do
        
    reads="${initDir}${f}_L001_R1_001.fastq.gz ${initDir}${f}_L002_R1_001.fastq.gz ${initDir}${f}_L003_R1_001.fastq.gz ${initDir}${f}_L004_R1_001.fastq.gz ${initDir}${f}_L001_R2_001.fastq.gz ${initDir}${f}_L002_R2_001.fastq.gz ${initDir}${f}_L003_R2_001.fastq.gz ${initDir}${f}_L004_R2_001.fastq.gz" 
    echo ${reads}
    
    sbatch -n ${nthread} --mem=16G -t 0-00:10:00 --job-name=${program1}${suffix}_${f} \
    -o ${baseDir}/LogErr/${program1}${suffix}_${f}.log -e ${baseDir}/LogErr/${program1}${suffix}_${f}.err \
    -p short --mail-user=${notifEm} --mail-type=ALL \
    --wrap="$program1 ${reads} --noextract -outdir ${baseDir}${program1}${suffix}/ > ${baseDir}${program1}${suffix}/fastQC_${f}.log;" 

done

echo 'end of fastqc script'

