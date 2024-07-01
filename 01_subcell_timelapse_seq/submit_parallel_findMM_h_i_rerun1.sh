#!/bin/bash
#SBATCH -n 1
#SBATCH -t 20       
#SBATCH --mem=10M                   
#SBATCH -p priority     

### UPDATED:	  10/23/2019 by Brendan Smalec
### LAST UPDATED: 11/15/2023 by Robert Ietswaart

### USE:			This script reruns findMismatches_complete_23_11_14_RI_variableGenome.R on files that did not finish in the first run
###					to output the read position, read nt, genome position, and genome nt
###					of each mismatch within a read. It also calculates the total number of
###					mismatches per read, the number of T>C mismatches within a read,
###					whether or not the read contains any T>C mismatches (binary call), and
###					finally the length of the read. 
###					Note that there has been several iterations of findMismatches_complete
###					and although we now convert reads to pseudo-bed format before running,
###					we used to run it with .sam files. Some of the comments and file 
###					naming may still refer to .sam format, however all code has been 
###					changed to run on .bed format. 

### REQUIREMENTS:	Processed alignments, split up into smaller files (10000 lines) with prefix tmp.


# HERE set the fasta path
fasta=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/K562_ensGRCh38_dm6_ercc_cat.fasta
#fasta=/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/NIH3T3_mm10_dm6_ercc_cat.fasta
#OLD fasta=/n/groups/churchman/bms36/genomes/hg38_dm6_ercc_grandslam/mtc_snp_masked/K562_ensGRCh38_dm6_ercc_cat.fasta
#OLD fasta=/n/groups/churchman/bms36/genomes/mm10_dm6_ercc_grandslam/mtc_snp_masked/mouseNIH3T3_mm10_dm6_ercc_cat.fasta

base_dir=`pwd`
Libs="T1_cyto T2_cyto T3_cyto T4_cyto T5_cyto" 
#Libs="T1_chr T2_chr T3_chr T4_chr T5_chr" 
#Libs="T1_poly T2_poly T3_poly T4_poly T5_poly" 
#Libs="U1_cyto U2_cyto U3_cyto U4_cyto U5_cyto"
#Libs="U1_chr U2_chr U3_chr U4_chr U5_chr"
#Libs="U1_poly U2_poly U3_poly U4_poly U5_poly"

#Libs=T1_cyto

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

	# LIST all tmp files to run through
	tmp_files=$(GLOBIGNORE="*.bed"; ls tmp*)
	unset GLOBIGNORE
	#tmp_files=tmp.ft

	# resubmit a job to run findMismatches_complete_23_11_14_RI_variableGenome.R on each tmp alignment file

	for x in $tmp_files
	do
	    if [ ! -f "${x}_MM.bed" ]; then
		OUT=${x}_mm.sh 		#THIS is the name of the script to be written. Include var x.
		echo "#!/bin/bash" > $OUT #NOTE: this overwrites any previous $OUT script to avoid appending of reruns
		echo "#SBATCH -n 1                              # Request one core" >> $OUT
		echo "#SBATCH -t 60" >> $OUT
		echo "#SBATCH -p short                           # Partition to run in" >> $OUT
		echo "#SBATCH --mem=500M" >> $OUT
		echo "" >> $OUT				
		echo "module load gcc/9.2.0" >> $OUT
		echo "module load bedtools/2.30.0" >> $OUT
		echo "module load R/4.3.1" >> $OUT
		echo "mkdir ${x}_analysis" >> $OUT
		echo "cp /n/groups/churchman/bms36/sandbox/TimelapseSeq/MismatchScripts_batchSubmission/findMismatches_complete_23_11_14_RI_variableGenome.R findMM_${x}.R" >> $OUT
		echo "./findMM_${x}.R ${x} ${current_dir} ${fasta}" >> $OUT
		echo "rm -r ${x}_analysis" >> $OUT
		echo "rm findMM_${x}.R" >> $OUT
		echo "" >> $OUT
		chmod u+x $OUT						# MUST change permission, otherwise will not submit
		sbatch $OUT							# Submit the script 
	    fi
	done
    done
done
echo "sbatch submit_parallel_findMM_h_i_rerun1.sh finished"

#Notes:

# once all tmp files have been analyzed, cat them together:
# cat tmp*MM.bed > reads.MM.bed
