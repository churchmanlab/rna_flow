#!/bin/bash
#Robert Ietswaart 
#sbatch job to proces GS top and bottom results into one output for fraction new RNA visualizations

date=20240120
org=m #h: human for K562, m: mouse for 3T3 
baseDir=/n/groups/churchman/ri23/bseq/Bayes${date}_3T3 #K562 for human, 3T3 for mouse

program=python
script=Timescale_Bayes_${date}_iii_prerun.py

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

jobname=${script}
    
sbatch -p priority -n ${nthread} --mem=1G -t 0-00:01:00 --job-name=${jobname} \
   -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
   --mail-user=${notifEm} --mail-type=NONE \
   --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org}"
