#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayesian rates calculation in batches.

date=20231213
org=h #h: human for K562, m: mouse for 3T3 
org_name=human #mouse
baseDir=/n/groups/churchman/ri23/bseq/BayesFactor4_${date}_K562 #K562 for human, 3T3 for mouse
input_iv=BayesFactor${date}_${org_name}_iv.log
N_BATCH=$(wc -l ${baseDir}/LogErr/${input_iv} | awk '{print $1}')
echo $N_BATCH

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
oflow=1e-100 #run 1/2: 1e-100, run 3: 1(?)
uflow=1 #run 1: 1e200, run 2: 1e100, run 3: 1(?) others?

for line in $(seq 2 $N_BATCH) #not line1: header
do    
    read -r prev_batch n_done start_id n_todo suffix_id <<<$(sed -n ${line}p ${baseDir}/LogErr/${input_iv})

    jobname=BayesFactor4_${date}_${org}_batch${suffix_id}

    if [ $n_todo -eq 1 ]
    then
        #retry the one that did not finish by itself in medium queue           
        sbatch -p short -n ${nthread} --mem=1G -t 0-12:00:00 --job-name=${jobname} \
            -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
            --mail-user=${notifEm} --mail-type=NONE \
            --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --oflow ${oflow} --uflow ${uflow} --start_id $start_id --n_genes $n_todo --suffix_id $suffix_id"

    else
        sbatch -p short -n ${nthread} --mem=1G -t 0-12:00:00 --job-name=${jobname} \
            -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
            --mail-user=${notifEm} --mail-type=NONE \
            --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --oflow ${oflow} --uflow ${uflow} --start_id $start_id --n_genes $n_todo --suffix_id $suffix_id"      
    fi
done

echo "end of Bayes_factor4_20231213_K562_iv.sh"
#short -n ${nthread} --mem=1G -t 0-12:00:00