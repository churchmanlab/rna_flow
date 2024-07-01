#!/bin/bash
#Robert Ietswaart 
#sbatch job to merge batches of Bayes Factor 4 compartment method

date=20240112
baseDir=/n/groups/churchman/ri23/bseq/BayesFactor20240112

program=python
script=Bayes_factor4_${date}_merge_with_BF3.py

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

jobname=${script}
    
sbatch -p priority -n ${nthread} --mem=1G -t 0-00:02:00 --job-name=${jobname} \
   -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
   --mail-user=${notifEm} --mail-type=NONE \
   --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script"
