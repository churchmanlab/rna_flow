#!/bin/bash

#SBATCH -c 1
#SBATCH -t 700
#SBATCH -p short
#SBATCH --mem=30G


### 20231121 Robert Ietswaart 

module load gcc/9.2.0
module load samtools/1.15.1
module load bedtools/2.30.0
module load R/4.3.1

base_dir=`pwd`

#Libs="T1_nuc T2_nuc T3_nuc T4_nuc T5_nuc T1_tot T2_tot T3_tot T4_tot T5_tot U1_nuc U2_nuc U3_nuc U4_nuc U5_nuc U1_tot U2_tot U3_tot U4_tot U5_tot"
Libs="T1_cyto T2_cyto T3_cyto T4_cyto T5_cyto T1_chr T2_chr T3_chr T4_chr T5_chr T1_poly T2_poly T3_poly T4_poly T5_poly" 
#Libs="U1_cyto U2_cyto U3_cyto U4_cyto U5_cyto U1_chr U2_chr U3_chr U4_chr U5_chr U1_poly U2_poly U3_poly U4_poly U5_poly"
#Libs='T1_nuc'

TURN_TYPES="_MM_temp_turnover _MM_temp_turnover_slow"

for lib in $Libs
do
    for suffix in $TURN_TYPES
    do
	cd ${base_dir}/${lib}/${lib}${suffix}
	current_dir=`pwd`

	# remove the previous scripts and any temp_analysis folders for each tmp file
	rm -f *_mm.sh
	rm -fr tmp.*_analysis

	# cat temp files
	cat tmp*MM.bed > reads.MM.bed

	# remove old files: not yet in case of need for rerun
	#rm -f submit_parallel_findMM.sh
	#rm -f test_submitMM.sh
	#rm slurm*
	#rm tmp*

	# sort according to read name 
	echo "${lib}${suffix}" 
	echo "Sorting according to read name"
	cp /n/groups/churchman/bms36/sandbox/TimelapseSeq/MismatchScripts_batchSubmission/sortReads_byName_23_11_19_RI.R srtR${lib}${suffix}.R
	./srtR${lib}${suffix}.R reads.MM.bed .
	rm srtR${lib}${suffix}.R
	echo "Done"


	# break into smaller files to run MM_per_fragment in parallel
	echo "Splitting up paired reads"
	split --lines=40000 reads.MM.bed_sort.bed rpf.
	echo "Done"

    done
done

echo "sbatch prerun_MM_per_fragment_h_ii.sh finished"
