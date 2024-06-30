#!/bin/bash
#Robert Ietswaart 
#sbatch job to run ML on subcellular rates for feature selection round 1.

date=20240122

program=python
script=ML_${date}_subcell_feat_select1.py

notifEm=robert_ietswaart@hms.harvard.edu
baseDir=/n/groups/churchman/ri23/bseq/ML${date}
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr
mkdir -p ${baseDir}/select_features

nthread=1 
     
N_FEAT_CLASSES=14 #see .py script

RATES="chr chr_release nuc nucexp cyto poly_entry whole_cell"

for r in $RATES
do
    for fc in $(seq 1 $N_FEAT_CLASSES)
    do     
        jobname=ML_${date}_subcell_select1_k_${r}_feat${fc}

        if [ $fc -ge 4 ] && [ $fc -le 6 ]; #kmers
        then
            mem_g=32G
            runtime="8:00:00"    
        else
            mem_g=16G
            runtime="1:00:00"
        fi
        
        sbatch -p short -n ${nthread} --mem=${mem_g} -t 0-${runtime} --job-name=${jobname} \
            -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
            --mail-user=${notifEm} --mail-type=ALL \
            --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --rate ${r} --feats ${fc}"

    done
done
