#!/bin/bash
#Robert Ietswaart 
#sbatch job to run ML on subcellular rates for feature selection round 2.

date=20230329

program=python
script=ML_${date}_subcell_feat_select2.py

notifEm=robert_ietswaart@hms.harvard.edu
baseDir=/n/groups/churchman/ri23/bseq/ML20230329
scDir=/n/groups/churchman/ri23/code/

module load gcc/6.2.0
module load python/3.7.4
module load java/jdk-1.8u112

cd ${scDir}
mkdir -p ${baseDir}/LogErr

nthread=1

RATES="chr chr_release nuc nucexp cyto poly_entry whole_cell"

for r in $RATES
do    
    jobname=ML_${date}_subcell_select2_k_${r}

    mem_g=16G
    runtime="4:00:00"    
        
    sbatch -p short -n ${nthread} --mem=${mem_g} -t 0-${runtime} --job-name=${jobname} \
    -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
    --mail-user=${notifEm} --mail-type=ALL \
    --wrap="source RNAdecayenv/bin/activate; cd bseq; $program $script --rate ${r}"

done
