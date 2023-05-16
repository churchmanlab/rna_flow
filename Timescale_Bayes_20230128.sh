#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayes factor calculation to distinguish residence and nucdeg models in batches.

date=20230128

program=python
script=Timescale_Bayes_${date}.py #for both mouse / human

notifEm=robert_ietswaart@hms.harvard.edu
baseDir=/n/groups/churchman/ri23/bseq/Bayes20230128 #for both mouse / human
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

org=m #h   
N_BATCH=0 #6051
#human: 6051
#mouse: 5861 
genes_w_rates_folder=Bayes20210615
#human: Bayes20210804_human
#mouse: Bayes20210615


for batch in $(seq 0 $N_BATCH)
do     
    jobname=Bayes${date}_timescales_${org}_batch${batch}
    
    sbatch -p short -n ${nthread} --mem=2G -t 0-12:00:00 --job-name=${jobname} \
        -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
        --mail-user=${notifEm} --mail-type=NONE \
        --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --start_id ${batch} --organism ${org} --gwr_folder ${genes_w_rates_folder}"

done
#medium -n ${nthread} --mem=1G -t 5-00:00:00
#short -n ${nthread} --mem=1G -t 0-12:00:00