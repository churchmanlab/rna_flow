#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayesian rates calculations for nucdeg rates .

date=20240120
org=h #h: human for K562, m: mouse for 3T3 
baseDir=/n/groups/churchman/ri23/bseq/Bayes${date}_K562 #K562 for human, 3T3 for mouse
N_BATCH=2497 #ceil((wc -l genes_w_nucdeg_reruns.tsv) / batch_size - 1) in $script)
#human: 2498  (wc -l genes_w_nucdeg_reruns.tsv) 
#mouse: tbd  (wc -l genes_w_nucdeg_reruns.tsv) 

program=python
script=Timescale_Bayes_${date}_nucdeg.py #for both mouse / human

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

#numerical over/underflow correction factors for the nucdeg rates posteriors
oflow=1e-100 #run 1: 1e-100, run 2: 1?,
uflow=1e250 #run 1/2: 1e100, run 3(?): 1e200, others?


for batch in $(seq 0 $N_BATCH)
do     
    jobname=Bayes${date}_rates_${org}_batch${batch}
    
    sbatch -p short -n ${nthread} --mem=1G -t 0-12:00:00 --job-name=${jobname} \
        -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
        --mail-user=${notifEm} --mail-type=NONE \
        --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --oflow ${oflow} --uflow ${uflow} --start_id ${batch}"

done
echo "end of Timescale_Bayes_20240120_K562_iii.sh"

#short -n ${nthread} --mem=1G -t 0-12:00:00
#priority -n ${nthread} --mem=1G -t 0-00:20:00
