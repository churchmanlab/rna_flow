#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayesian rates calculation in batches.

date=20231213
org=m #h: human for K562, m: mouse for 3T3 
baseDir=/n/groups/churchman/ri23/bseq/BayesFactor4_${date}_3T3 #K562 for human, 3T3 for mouse
N_BATCH=1798 #floor((wc -l genes_w_rates.csv) / batch_size in Timescale_Bayes_${date}.py script)
#human:  1816 (wc -l genes_w_rates.csv) = 18165, current batch_size=10, or revert to 2.
#mouse: 1798 (wc -l genes_w_rates.csv) = 17985

program=python
script=Bayes_factor4_${date}.py #for both mouse / human

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

oflow=1e-100 #run 1/2: 1e-100, run 3: 1(?)
uflow=1e200 #run 1: 1e200, run 2: 1e100, run 3: 1(?) others?


for batch in $(seq 0 $N_BATCH)
do     
    jobname=BayesFactor4_${date}_${org}_batch${batch}
    
    sbatch -p short -n ${nthread} --mem=1G -t 0-12:00:00 --job-name=${jobname} \
        -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
        --mail-user=${notifEm} --mail-type=NONE \
        --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --oflow ${oflow} --uflow ${uflow} --start_id ${batch}"

done
echo "end of Bayes_factor4_20231213_3T3_i.sh"
#short -n ${nthread} --mem=1G -t 0-12:00:00