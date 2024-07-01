#!/bin/bash
#Robert Ietswaart 
#sbatch job to estimate TC conversion rates (pc, perr1, perr2) for each human K562 sample

date=20231125

program=python
script=TCconversion_from_background_${date}_K562.py

notifEm=robert_ietswaart@hms.harvard.edu
baseDir=/n/groups/churchman/ri23/bseq/TC20231125_K562
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

   
jobname=TCconversion_from_background_${date}_K562
    
sbatch -p short -n ${nthread} --mem=16G -t 0-12:00:00 --job-name=${jobname} \
    -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
    --mail-user=${notifEm} --mail-type=ALL \
    --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script"