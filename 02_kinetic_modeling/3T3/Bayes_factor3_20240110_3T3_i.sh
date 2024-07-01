#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayes Factor calculation using 3 compartments in batches.

date=20240110
org=m #h: human for K562, m: mouse for 3T3 
baseDir=/n/groups/churchman/ri23/bseq/BayesFactor3_${date}_3T3 #K562 for human, 3T3 for mouse
N_BATCH=719 #floor((wc -l genes_w_rates.csv) / batch_size in Timescale_Bayes_${date}.py script)
#human:  726 (wc -l genes_w_rates.csv) = 18165, current batch_size=25.
#mouse: 719 (wc -l genes_w_rates.csv) = 17985

program=python
script=Bayes_factor3_${date}.py #for both mouse / human

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

for batch in $(seq 0 $N_BATCH)
do     
    jobname=BayesFactor3_${date}_${org}_batch${batch}
    
    sbatch -p short -n ${nthread} --mem=1G -t 0-12:00:00 --job-name=${jobname} \
        -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
        --mail-user=${notifEm} --mail-type=NONE \
        --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --start_id ${batch}"

done
echo "end of Bayes_factor3_20240110_3T3_i.sh"
#short -n ${nthread} --mem=1G -t 0-12:00:00