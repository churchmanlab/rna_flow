#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayesian analysis in batches

date=20210804 

program=python
script=Timescale_fit_20210712_human.py

notifEm=robert_ietswaart@hms.harvard.edu
baseDir=/n/groups/churchman/ri23/bseq/MLE20210712_human
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 
     
jobname=Timescale_fit_${date}_human

sbatch -p short -n ${nthread} --mem=16G -t 0-12:00:00 --job-name=${jobname} \
    -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
    --mail-user=${notifEm} --mail-type=ALL \
    --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script > ${baseDir}/LogErr/${jobname}.out"


#short -n ${nthread} --mem=16G -t 0-00:10:00        #test run
#short -n ${nthread} --mem=16G -t 0-12:00:00    #full batch runs