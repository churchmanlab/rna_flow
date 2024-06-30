#!/bin/bash
#Robert Ietswaart 
#sbatch job to align Nuc Wash fastq files with STAR 

date=20231005
scDir=/n/groups/churchman/ri23/code/bseq/
cd ${scDir}
suffix="" #default: "" #also tested: _barcode-mismatches0
parameterFiles=NucWash_CG20230908${suffix}_parameters.in #
echo $parameterFiles 

baseDir=`grep "Base Directory"            $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
samples=`grep "Sample Names"           $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
indexDir=`grep "Index Directory"            $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
initDir=`grep "Initial Files Directory"     $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
notifEm=`grep "Notification Email"       $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`

program1=STAR
program2=samtools

module load gcc/6.2.0
module load star/2.7.0a
module load ${program2}/1.15.1

pFile=${scDir}/NucWash_STAR_parameters.in
nthread=4
c=_cutadapt #if no cutadapt adapter trimming (not used): ""

FILES=`echo -e $samples`
for f in $FILES
do
    echo "Doing file "$f
    outDir=${baseDir}STAR/${f}
    mkdir -p ${outDir}/LogErr
    cd ${outDir}
    reads1=${initDir}${f}_L001_R1_001${c}.fastq.gz,${initDir}${f}_L002_R1_001${c}.fastq.gz,${initDir}${f}_L003_R1_001${c}.fastq.gz,${initDir}${f}_L004_R1_001${c}.fastq.gz
    reads2=${initDir}${f}_L001_R2_001${c}.fastq.gz,${initDir}${f}_L002_R2_001${c}.fastq.gz,${initDir}${f}_L003_R2_001${c}.fastq.gz,${initDir}${f}_L004_R2_001${c}.fastq.gz   

    #Align with STAR
    #Sort and index resulting bam files with samtools
    sbatch -n ${nthread} --mem-per-cpu=16G -t 0-01:00:00 --job-name=${program1}_${f}${c} \
    -o ${outDir}/LogErr/${program1}_${f}${c}.log -e ${outDir}/LogErr/${program1}_${f}${c}.err \
    -p short --mail-user=${notifEm} --mail-type=ALL \
    --wrap="$program1 --genomeDir ${indexDir} --readFilesIn ${reads1} ${reads2} \
    --runThreadN ${nthread} --parametersFiles ${pFile} --outFileNamePrefix ${f}_ ; \
    $program2 sort -@ ${nthread} -o ${f}_sorted_indexed.bam ${f}_Aligned.sortedByCoord.out.bam ; \
    $program2 index -@ ${nthread} ${f}_sorted_indexed.bam ; " 
done

echo 'end of bulk RNAseq_align sh script'

#sbatch -n ${nthread} --mem-per-cpu=16G -t 0-01:00:00
#2B_NW_NRRD_S18 TIMEOUT -> NucWash_CG20230908${suffix}_parameters_rerun.in #: sbatch -n ${nthread} --mem-per-cpu=16G -t 0-01:30:00 
