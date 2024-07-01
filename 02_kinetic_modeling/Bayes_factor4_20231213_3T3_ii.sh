#!/bin/bash
#Robert Ietswaart 
#sbatch job to run Bayesian rates calculation in batches.

date=20231213
org=m #h: human for K562, m: mouse for 3T3 
org_name=mouse #human #
baseDir=/n/groups/churchman/ri23/bseq/BayesFactor4_${date}_3T3 #K562 for human, 3T3 for mouse
N_BATCH=1798 #floor((wc -l genes_w_rates.csv) / batch_size in Timescale_Bayes_${date}.py script)
#human:  1816 (wc -l genes_w_rates.csv) = 18165, current batch_size=10, or revert to 2.
#mouse: 1798 (wc -l genes_w_rates.csv) = 17985

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
batch_size=10 #has to match value in .py script: 10, or revert to 2.
oflow=1e-100 #run 1/2: 1e-100, run 3: 1(?)
uflow=1e200 #run 1: 1e200, run 2: 1e100, run 3: 1(?) others?

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
            
            jobname=BayesFactor4_${date}_${org}_batch${batch}_ii_first
            sbatch -p medium -n ${nthread} --mem=1G -t 5-00:00:00 --job-name=${jobname} \
                -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
                --mail-user=${notifEm} --mail-type=NONE \
                --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --oflow ${oflow} --uflow ${uflow} --start_id $start_id --n_genes $n_todo --suffix_id $suffix_id"

            echo $batch, $n_done, $start_id, $n_todo, $suffix_id >> ${baseDir}/LogErr/BayesFactor20231213_${org_name}_ii.log
            suffix_id=$(($suffix_id+1))

            #then set the rest done
            n_done=1
        fi #n_done==0
               
        start_id=$(($batch*$batch_size+$n_done))
        n_todo=$(($batch_size-$n_done))
       
        jobname=BayesFactor4_${date}_${org}_batch${batch}_ii_rest
        sbatch -p short -n ${nthread} --mem=1G -t 0-12:00:00 --job-name=${jobname} \
            -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
            --mail-user=${notifEm} --mail-type=NONE \
            --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --organism ${org} --oflow ${oflow} --uflow ${uflow} --start_id $start_id --n_genes $n_todo --suffix_id $suffix_id"
           
        echo $batch, $n_done, $start_id, $n_todo, $suffix_id >> ${baseDir}/LogErr/BayesFactor20231213_${org_name}_ii.log
        suffix_id=$(($suffix_id+1))
       
    fi #rerun required: n_done < batch_size
    
done





echo "end of Bayes_factor4_20231213_3T3_ii.sh"
#short -n ${nthread} --mem=1G -t 0-12:00:00