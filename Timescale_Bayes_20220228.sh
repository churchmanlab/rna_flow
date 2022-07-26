#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayesian (residence and nucdeg) models in batches with improvements.

date=20220228

program=python
script=Timescale_Bayes_${date}_human.py
#mouse: Timescale_Bayes_${date}.py 
#human: Timescale_Bayes_${date}_human.py

notifEm=robert_ietswaart@hms.harvard.edu
baseDir=/n/groups/churchman/ri23/bseq/Bayes${date}_human #mouse: Bayes20220222 #
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 
     
N_BATCH=181
#mouse: 175 
#human: 181 

for batch in $(seq 0 $N_BATCH)
do     
    jobname=Bayes${date}_human_batch${batch}
    #mouse: Bayes${date}_batch${batch} 
    #human: Bayes${date}_human_batch${batch}
    
    sbatch -p short -n ${nthread} --mem=16G -t 0-12:00:00 --job-name=${jobname} \
        -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
        --mail-user=${notifEm} --mail-type=ALL \
        --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --start_id ${batch}"

done

#short -n ${nthread} --mem=16G -t 0-12:00:00    #full batch runs