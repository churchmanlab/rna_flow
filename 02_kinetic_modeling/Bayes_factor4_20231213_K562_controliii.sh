#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayes Factor calculation in batches.

date=20231213
org=h #h: human for K562, m: mouse for 3T3 
baseDir=/n/groups/churchman/ri23/bseq/BayesFactor4_${date}_K562_oflow1 #K562 for human, 3T3 for mouse
start_id=6651 #0-based, 6652 1-based #grep -n ENSG00000135486 genes_w_rates.csv 
N_BATCH=10002 
n_todo=1


program=python
script=Bayes_factor4_${date}_rerun.py #for both mouse / human

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

oflow=1 #control run 1: to see if the results change compared to 1e-100
uflow=1e300 #run 1: 1e200, run 2: 1e100, run 3: 1(?) others?


for suffix_id in $(seq $N_BATCH $N_BATCH)
do     
    jobname=BayesFactor4_${date}_${org}_batch${suffix_id}
    
    sbatch -p short -n ${nthread} --mem=1G -t 0-12:00:00 --job-name=${jobname} \
        -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
        --mail-user=${notifEm} --mail-type=NONE \
        --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --oflow ${oflow} --uflow ${uflow} --start_id $start_id --n_genes $n_todo --suffix_id $suffix_id"

done
echo "end of Bayes_factor4_20231213_K562_controliii.sh"
#short -n ${nthread} --mem=1G -t 0-12:00:00