#!/bin/bash
#Robert Ietswaart 
#sbatch job to merge batches of Bayes Factor 4 compartment method

date=20240110
org=m #h: human for K562, m: mouse for 3T3 
baseDir=/n/groups/churchman/ri23/bseq/BayesFactor3_${date}_3T3 #K562 for human, 3T3 for mouse
nbatch=1173 
#h: 1237+1=1238 
#m: 1172+1=1173 

program=python
script=Bayes_factor3_${date}_merge.py

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

jobname=${script}
    
sbatch -p priority -n ${nthread} --mem=1G -t 0-00:05:00 --job-name=${jobname} \
   -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
   --mail-user=${notifEm} --mail-type=NONE \
   --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --nbatch $nbatch"
