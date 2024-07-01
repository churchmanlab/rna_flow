#!/bin/bash
#SBATCH -c 1
#SBATCH -t 60
#SBATCH --mem=100M
#SBATCH -p priority

### 20231119 Robert Ietswaart
### Batch submit version of run_MM_per_fragment.sh (by Brendan Smalec)

### USE: 			This script will concatenate the mismatches from paired reads into 
###					one line (total mismatches per fragment) while also removing 
###					mismatches that are duplicated on both reads (aka, they are within
###					the overlapping region). 

###					The MM_per_fragment script keeps the information for read #1 of the
###					pair (coordinates, sequence, etc) while updating the following fields
###					to reflect the mismatch information for both read #1 and read #2:
###					mismatch read nt, mismatch genome, position, mismatch genome nt, 
###					total # of mismatches, # of T>C mismatches, and whether or not the 
###					fragment contains any T>C mismatches (binary call). However this 
###					does NOT update the mismatch read position or the length of the read,
###					so if these parameters will be used for subsequent analysis, the 
###					script will need to be adapted. 

### REQUIREMENTS:	                rpf.*.MM.bed 
###					MM_per_fragment.R 

# load required modules
module load gcc/9.2.0
module load samtools/1.15.1
module load bedtools/2.30.0
module load R/4.3.1

base_dir=`pwd`

#Libs="T1_nuc T2_nuc T3_nuc T4_nuc T5_nuc T1_tot T2_tot T3_tot T4_tot T5_tot U1_nuc U2_nuc U3_nuc U4_nuc U5_nuc U1_tot U2_tot U3_tot U4_tot U5_tot"
Libs="T1_cyto T2_cyto T3_cyto T4_cyto T5_cyto T1_chr T2_chr T3_chr T4_chr T5_chr T1_poly T2_poly T3_poly T4_poly T5_poly U1_cyto U2_cyto U3_cyto U4_cyto U5_cyto U1_chr U2_chr U3_chr U4_chr U5_chr U1_poly U2_poly U3_poly U4_poly U5_poly"


TURN_TYPES="_MM_temp_turnover _MM_temp_turnover_slow"

for lib in $Libs
do
    for suffix in $TURN_TYPES
    do
	cd ${base_dir}/${lib}/${lib}${suffix}
	current_dir=`pwd`

	# remove the previous scripts and any temp_analysis folders for each rpf file
	rm -f *_MMperF.sh
	rm -fr rpf.*_analysis

	# LIST all tmp files to run through
	tmp_files=$(GLOBIGNORE="*.bed"; ls rpf*)
	unset GLOBIGNORE

	# submit all tmp files for analysis in parallel
	for x in $tmp_files
	do
	
		OUT=${x}_MMperF.sh 		#THIS is the name of the script to be written. Include var x.
		echo "#!/bin/bash" > $OUT       #NOTE: this overwrites any previous $OUT script to avoid appending of reruns
		echo "#SBATCH -n 1                              # Request one core" >> $OUT
		echo "#SBATCH -t 30" >> $OUT   #120, TEMP trial of 30 for all
		echo "#SBATCH -p short                           # Partition to run in" >> $OUT
		echo "#SBATCH --mem=500M" >> $OUT
		echo "" >> $OUT				
		echo "module load gcc/9.2.0" >> $OUT
		echo "module load R/4.3.1" >> $OUT
		echo "mkdir ${x}_analysis" >> $OUT
		echo "cd ${x}_analysis" >> $OUT
		echo "" >> $OUT
		echo "cp /n/groups/churchman/bms36/sandbox/TimelapseSeq/MismatchScripts_batchSubmission/MM_per_fragment_23-11-19_RI_batchSubmit.R MMpf_${x}.R" >> $OUT
		echo "./MMpf_${x}.R ${x} ${current_dir} ${x}" >> $OUT
		echo "mv ${x}_fragments.bed ../" >> $OUT
		echo "cd .." >> $OUT
		echo "rm -r ${x}_analysis" >> $OUT
		echo "" >> $OUT
		chmod u+x $OUT						# MUST change permission, otherwise will not submit
		sbatch $OUT							# Submit the script 
	
	done
    done
done

echo "sbatch submit_run_MM_per_fragment_h_ii.sh finished"

# when this is done running, test that all temporary read per fragment (rpf*) files have been analyzed by running:
# ./test_MM_per_F.sh
# like with test_submitMM.sh

# once all of these have run, then run in a next script:
# cat rpf*fragments.bed > reads_MM.fragments.bed


