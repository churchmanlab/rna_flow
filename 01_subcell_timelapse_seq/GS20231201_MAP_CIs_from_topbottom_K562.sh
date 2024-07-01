#!/bin/bash
#Robert Ietswaart 
#sbatch job to proces GS top and bottom results into one output for fraction new RNA visualizations

date=20231201
org=h #h: human for K562, m: mouse for 3T3 
baseDir=/n/groups/churchman/ri23/bseq/GS${date}_K562 #3T3

program=python
script=GS${date}_MAP_CIs_from_topbottom.py

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

jobname=${script}
    
sbatch -p short -n ${nthread} --mem=2G -t 0-06:00:00 --job-name=${jobname} \
   -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
   --mail-user=${notifEm} --mail-type=NONE \
   --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org}"