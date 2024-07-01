#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayesian analysis in batches

date=20240118 
org=h #h: human for K562, m: mouse for 3T3 

program=python
script=Timescale_LSE_${date}.py

notifEm=robert_ietswaart@hms.harvard.edu
baseDir=/n/groups/churchman/ri23/bseq/LSE20240118
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 
     
jobname=Timescale_LSE_${date}_${org}

sbatch -p priority -n ${nthread} --mem=2G -t 0-12:00:00 --job-name=${jobname} \
    -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
    --mail-user=${notifEm} --mail-type=NONE \
    --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org}"


#priority -n ${nthread} --mem=1G -t 0-00:01:00        #test run
#priority -n ${nthread} --mem=2G -t 0-12:00:00    #full batch runs