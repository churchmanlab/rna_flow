#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayes factor calculation to distinguish residence and nucdeg models in batches.

date=20220315

program=python
script=Bayes_factor_${date}.py #for both mouse / human

notifEm=robert_ietswaart@hms.harvard.edu
baseDir=/n/groups/churchman/ri23/bseq/BayesFactor20220307 #for both mouse / human
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

org=h #m 
N_BATCH=726
#mouse: 703 
#human: 726 
genes_w_rates_folder=Bayes20210804_human
#mouse: Bayes20210615
#human: Bayes20210804_human


for batch in $(seq 0 $N_BATCH)
do     
    jobname=BayesFactor${date}_${org}_batch${batch}
    
    sbatch -p short -n ${nthread} --mem=16G -t 0-12:00:00 --job-name=${jobname} \
        -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
        --mail-user=${notifEm} --mail-type=ALL \
        --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --start_id ${batch} --organism ${org} --gwr_folder ${genes_w_rates_folder}"

done
#short -n ${nthread} --mem=16G -t 0-12:00:00    #full batch runs