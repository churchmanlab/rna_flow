#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayesian rates calculation in batches.

date=20240119
org=h #h: human for K562, m: mouse for 3T3 
baseDir=/n/groups/churchman/ri23/bseq/Bayes${date}_K562 #K562 for human, 3T3 for mouse
N_BATCH=1816 #floor((wc -l genes_w_rates.csv) / batch_size in Timescale_Bayes_${date}.py script)
#human:  1816 (wc -l genes_w_rates.csv) = 18165
#mouse: 1798 (wc -l genes_w_rates.csv) = 17985

program=python
script=Timescale_Bayes_${date}.py #for both mouse / human

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

#numerical over/underflow correction factors for the nucdeg rates posteriors
oflow=1e-100 #run 1: 1e-100, run 2: 1e-200/1e-100, run 3: 1e-100
uflow=1e100 #run 1/2: 1e100, run 3: 1e250


for batch in $(seq 0 $N_BATCH)
do     
    jobname=Bayes${date}_rates_${org}_batch${batch}
    
    sbatch -p short -n ${nthread} --mem=1G -t 0-12:00:00 --job-name=${jobname} \
        -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
        --mail-user=${notifEm} --mail-type=NONE \
        --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --oflow ${oflow} --uflow ${uflow} --start_id ${batch}"

done
echo "end of Timescale_Bayes_${date}_K562_i.sh"
#short -n ${nthread} --mem=1G -t 0-12:00:00
