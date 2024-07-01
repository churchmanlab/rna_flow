#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayesian rates calculation in batches.

date=20240110
org=m #h: human for K562, m: mouse for 3T3 
org_name=mouse #human
baseDir=/n/groups/churchman/ri23/bseq/BayesFactor3_${date}_3T3 #K562 for human, 3T3 for mouse
N_BATCH=719 #floor((wc -l genes_w_rates.csv) / batch_size in Timescale_Bayes_${date}.py script)
#human:  726 (wc -l genes_w_rates.csv) = 18165, current batch_size=25, or revert to 10.
#mouse: 719 (wc -l genes_w_rates.csv) = 17985

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
batch_size=25 #has to match value in .py script: 10, or revert to 2.


suffix_id=$(($N_BATCH+1))


for batch in $(seq 0 $N_BATCH)
do     
    #complete file has n_done=batch_size+1 lines: 1 header row + batch_size genes 
    n_done=$(wc -l ${baseDir}/Bayes_factor_${date}_${org_name}_batch${batch}.tsv | awk '{print $1 - 1}')

    if [ $n_done -lt $batch_size ]
    then
        
        if [ $n_done -eq 0 ]
        then
        #retry the one that did not finish by itself in medium queue (should not happen often)
            
            start_id=$(($batch*$batch_size+$n_done))
            n_todo=1
            
            jobname=BayesFactor3_${date}_${org}_batch${suffix_id}
            sbatch -p medium -n ${nthread} --mem=1G -t 5-00:00:00 --job-name=${jobname} \
                -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
                --mail-user=${notifEm} --mail-type=NONE \
                --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --start_id $start_id --n_genes $n_todo --suffix_id $suffix_id"

            echo $batch, $n_done, $start_id, $n_todo, $suffix_id >> ${baseDir}/LogErr/BayesFactor${date}_${org_name}_ii.log
            suffix_id=$(($suffix_id+1))

            #then set the rest done
            n_done=1
        fi #n_done==0
               
        start_id=$(($batch*$batch_size+$n_done))
        n_todo=$(($batch_size-$n_done))
       
        jobname=BayesFactor3_${date}_${org}_batch${suffix_id}
        sbatch -p short -n ${nthread} --mem=1G -t 0-12:00:00 --job-name=${jobname} \
            -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
            --mail-user=${notifEm} --mail-type=NONE \
            --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --start_id $start_id --n_genes $n_todo --suffix_id $suffix_id"
           
        echo $batch, $n_done, $start_id, $n_todo, $suffix_id >> ${baseDir}/LogErr/BayesFactor${date}_${org_name}_ii.log
        suffix_id=$(($suffix_id+1))
       
    fi #rerun required: n_done < batch_size
    
done


echo "end of Bayes_factor3_${date}_3T3_ii.sh"
#short -n ${nthread} --mem=1G -t 0-12:00:00