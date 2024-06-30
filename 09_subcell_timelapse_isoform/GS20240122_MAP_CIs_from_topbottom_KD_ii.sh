#!/bin/bash
#Robert Ietswaart 
#sbatch job to proces GS top and bottom results into one output for fraction new RNA visualizations

date=20240122
suf=ii
baseDir=/n/groups/churchman/ri23/bseq/GS${date}_KD_${suf}

program=python
script=GS${date}_MAP_CIs_from_topbottom_KD_${suf}.py

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

read_types='unspliced_junctions introns not_introns'


for rt in $read_types
do
    jobname=${script}_${rt}

    sbatch -p short -n ${nthread} --mem=1G -t 0-01:00:00 --job-name=${jobname} \
       -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
       --mail-user=${notifEm} --mail-type=NONE \
       --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --read_type ${rt}"
  
done
#-p short -n ${nthread} --mem=1G -t 0-03:00:00
