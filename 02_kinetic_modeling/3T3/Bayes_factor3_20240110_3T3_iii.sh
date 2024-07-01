#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayesian rates calculation in batches.

date=20240110
org=m #h: human for K562, m: mouse for 3T3 
org_name=mouse #human
baseDir=/n/groups/churchman/ri23/bseq/BayesFactor3_${date}_3T3 #K562 for human, 3T3 for mouse
input_iii=BayesFactor${date}_${org_name}_iii.log
N_BATCH=$(wc -l ${baseDir}/LogErr/${input_iii} | awk '{print $1}')
echo $N_BATCH

program=python
script=Bayes_factor3_${date}_rerun.py #for both mouse / human

notifEm=robert_ietswaart@hms.harvard.edu
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1 

for line in $(seq 2 $N_BATCH) #not line1: header
do    
    read -r prev_batch n_done start_id n_todo suffix_id <<<$(sed -n ${line}p ${baseDir}/LogErr/${input_iii})

    jobname=BayesFactor3_${date}_${org}_batch${suffix_id}

    if [ $n_todo -eq 1 ]
    then
        #retry the one that did not finish by itself in medium queue           
        sbatch -p medium -n ${nthread} --mem=1G -t 5-00:00:00 --job-name=${jobname} \
            -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
            --mail-user=${notifEm} --mail-type=NONE \
            --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --start_id $start_id --n_genes $n_todo --suffix_id $suffix_id"

    else
        sbatch -p short -n ${nthread} --mem=1G -t 0-12:00:00 --job-name=${jobname} \
            -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
            --mail-user=${notifEm} --mail-type=NONE \
            --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --start_id $start_id --n_genes $n_todo --suffix_id $suffix_id"      
    fi
done

echo "end of Bayes_factor3_20240110_3T3_iii.sh"
#short -n ${nthread} --mem=1G -t 0-12:00:00