#!/bin/bash
#Robert Ietswaart 
#sbatch job to run bcl2fastq on Illumina NEXTseq raw output files 
#more info: https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html
#more info: bcl2fastq2 software User Guide: bcl2fastq2-v2-20-software-guide-15051736-03.pdf

date=20230912

program=bcl2fastq

notifEm=robert_ietswaart@hms.harvard.edu
baseDir=/n/groups/churchman/ri23/bseq/CG20230908/
seqDirs=230908_NB500910_0732_AHTHMHBGXT #YYMMDD_InstrumentID_RunID_FlowcellID
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load ${program}/2.20.0.422

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=4 

for batch in $seqDirs
do     
    # bcl2fastq default value: --barcode-mismatches 1
    jobname=${program}_${date}_${batch}
    echo ${program}_${date}_${batch}

    sbatch -p short -n ${nthread} --mem=32G -t 0-12:00:00 --job-name=${jobname} \
        -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
        --mail-user=${notifEm} --mail-type=ALL \
        --wrap="cd bseq; $program --runfolder-dir ${baseDir}$batch --output-dir ${baseDir}fastq --sample-sheet ${baseDir}CG20230908_NEXTseq_SampleSheet.csv -r $nthread -w $nthread -p $nthread --barcode-mismatches 1"
        
    #control: stringent, no mismatch allowed. 
    jobname=${program}_${date}_${batch}_barcode-mismatches0
    echo ${program}_${date}_${batch}

    sbatch -p short -n ${nthread} --mem=32G -t 0-12:00:00 --job-name=${jobname} \
        -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
        --mail-user=${notifEm} --mail-type=ALL \
        --wrap="cd bseq; $program --runfolder-dir ${baseDir}$batch --output-dir ${baseDir}fastq_barcode-mismatches0 --sample-sheet ${baseDir}CG20230908_NEXTseq_SampleSheet.csv -r $nthread -w $nthread -p $nthread --barcode-mismatches 0"

done

echo "done submissions" 
