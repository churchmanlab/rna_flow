# RNA flow

Code to accompany publication:  
Ietswaart, R., Smalec B.M., Xu A., et al, (2024).  
Genome-wide quantification of RNA flow across subcellular compartments reveals determinants of the mammalian transcript life cycle.  
Mol. Cell. 84, 1-20  
https://doi.org/10.1016/j.molcel.2024.06.008  




# 01_subcell_timelapse_seq

### Trim reads with Cutadapt
cutadapt.sh

### STAR alignment
make_STAR_dirs.sh  
run_align_human_ii.sh  
which calls: align_human.sh  

### GRAND-SLAM (with default parameters):
bamToCIT_human_ii.sh  
run_grandslam_human_ii_12combined.sh  

### Pc estimation
1. 
submit_process_n_k_turnover_WT_2023_h.sh  
which calls:  
process_n_k_turnover_WT_2023_h.sh  
which calls:  
top1000genes_turnover.R  
modifyBed.R  
MismatchScripts_batchSubmission/sortReads_byName_noSingles.R  
bottom500genes_turnover.R  

2. 
submit_parallel_findMM_h_i_rerun1.sh  
which calls: findMismatches_complete_23_11_14_RI_variableGenome.R  

3. 
prerun_MM_per_fragment_h_i.sh  
which calls: sortReads_byName_23_11_19_RI.R  

4. 
submit_run_MM_per_fragment_[m|h]_[i|ii|iii|iv].sh  
which calls: MM_per_fragment_23-11-19_RI_batchSubmit.R  

5. In the folders *MM_temp_turnover*  
run_shorten_fragments.sh  

6. Generate the T>C mismatch frequency matrices.  
run_MM_freq.sh  
which calls:  
MM_[rep]_[frac].R  
MM_[rep]_[frac]_slow.R  

7. Estimate the pc using the binomial mixture model  
TCconversion_from_background_20231125_K562.sh  
which calls: TCconversion_from_background_20231125_K562.py  

### Rerun GRAND-SLAM with pc / pe estimates from binomial mixture model for top and bottom.  
copy the GS_run1 folder to new folders top and bottom
remove the previous GS output files: [sample].tsv and [sample].ext.tsv
in top/bottom folders: replace the pc values with our respective top/bottom estimates into *.rates.tsv file in the rows single_new and double_new
then
For the -4sU timepoint (sample *1): replace the pc values with the 15min timepoint for consistency (it does not matter what it is though, since we do not analyse that timepoint).
remove the GS output files  
run_grandslam.sh  

Merge GS top and bottom results (used for visualizations in MS)  
GS20231201_MAP_CIs_from_topbottom_K562.sh  
which calls: GS20231201_MAP_CIs_from_topbottom.py  


# 02_kinetic_modeling
0. Generate a list of genes with GS data to analyze  
Timescale_Bayes_20231201_prerun.sh  
which calls: Timescale_Bayes_20231201_prerun.py  

1. Batch run i  
Timescale_Bayes_20240119_K562_i.sh  
which calls: Timescale_Bayes_20240119.py  

2. Batch run ii  
Timescale_Bayes_20240119_K562_ii.sh  
which calls:  
Timescale_Bayes_20240119_rerun.py  
Timescale_Bayes_20240119_rerun_no_nucdeg.py  

3.0. Merge batch rate files  
Timescale_Bayes_20240119_merge.sh  
which calls: Timescale_Bayes_20240191_merge.py  

3.1. Identify the genes that need a rerun for nucdeg rates  
Timescale_Bayes_20240120_iii_prerun.sh  
which calls: Timescale_Bayes_20240120_iii_prerun.py  

3.2 Run iii for nucdeg rates  
Timescale_Bayes_20240120_K562_iii.sh  
which calls: Timescale_Bayes_20240120_nucdeg.py  

4. Merge missing nucdeg rates into rate file  
Timescale_Bayes_20240120_merge_nucdeg.sh  
which calls: Timescale_Bayes_20240120_merge_nucdeg.py  

### Bayes Factor calculations
Four compartment Bayes factor calculation
Bayes_factor4_20231213_K562_i.sh  
which calls: Bayes_factor4_20231213.py  

Bayes_factor4_20231213_K562_ii.sh  
which calls: Bayes_factor4_20231213_rerun.py  

Bayes_factor4_20231213_iii_prerun.sh  
which calls: Bayes_factor4_20231213_iii_prerun.py  

Bayes_factor4_20231213_K562_iii.sh  
which calls: Bayes_factor4_20231213_iii_prerun.py  

Bayes_factor4_20231213_iv_prerun.sh  
which calls: Bayes_factor4_20231213_iv_prerun.py  

Bayes_factor4_20231213_K562_iv.sh  
which calls: Bayes_factor4_20231213_rerun.py  

Merge BF4 batches:  
Bayes_factor4_20231213_merge.sh  
which calls: Bayes_factor4_20231213_merge.py  

Three compartment Bayes factor calculation
Bayes_factor3_20240110_K562_i.sh  
which calls: Bayes_factor3_20240110.py  

Bayes_factor3_20240110_K562_ii.sh  
which calls: Bayes_factor3_20240110_rerun.py  

Bayes_factor3_20240110_iii_prerun.sh  
which calls: Bayes_factor3_20240110_iii_prerun.py  

Bayes_factor3_20240110_K562_iii.sh  
which calls: Bayes_factor3_20240110_rerun.py  

Merge BF3 batches:  
Bayes_factor3_20240110_merge.sh  
which calls: Bayes_factor3_20240110_merge.py  

Merge BF3 and BF4 results:  
Bayes_factor4_20240112_merge_with_BF3.sh  
which calls: Bayes_factor4_20240112_merge_with_BF3.py  


# 03_nuclear_wash_subcell_rna_seq
Convert raw Illumina seq output (.bcl) to fastq  
bcl2fastq_20230912.sh

fastQC (optional)  
NucWash_fastqc_20230920.sh

Hardclip Illumina adapter sequences  
NucWash_cutadapt_20231002.sh

Generate and concatenate human and yeast (spike-in) genome fasta and gtf / build STAR index  
NucWash_STAR_index_20230915.sh

run STAR with default paired-end RNAseq parameters  
NucWash_STAR_alignment_20231005.sh  
which calls:  
NucWash_CG20230908_parameters.in  
NucWash_STAR_parameters.in  

Get gene counts from alignments using featurecounts  
NucWash_featureCounts_20231005.sh  

Confirm that yeast spike in transcripts are in alignments (optional) 
NucWash_inspect_bam_20231005.sh  

Perform DE and PCA with spike in  
NucWash_DEseq2_20240203.Rmd

Boxplots nuclear half lives vs Nuc Wash RNA levels  
Nuc_cyto_ratios_NucWash_vs_rates_20240129.Rmd  

# 04_subcell_nanostring
Erikstrings_fit_20220328.ipynb

# 05_gi_map
Alignment:  
screen_processing_scripts_python3.sh

Analysis:  
AXSNF03_analysis_cleaned.Rmd  

# 06_ribosome_profiling
Alignment:  
ribosome-profiling-preprocess-alignment_all_scripts.sh  

Analysis:  
AXSNF03_ribosome-profiling-analysis-cleaned.Rmd

# 07_DDX3X_subcell_rna_seq
Alignment:  
fractionation_alignment_all_scripts.sh

Analysis:  
AXSNF10-analysis-cleaned.Rmd

# 08_subcell_nanopore_seq
RNA_flow_alignment_directRNA_ONT.sh

# 09_subcell_timelapse_isoform

Script to generate a modified gtf for use with GRAND-SLAM:  
GS20240122_Get_gene_level_gtf.sh

Script to generate a GRAND-SLAM index (.oml) file  
GS20240122_Get_gedi_oml.sh

Scripts to filter all KDii samples and run GRAND-SLAM  
GS20240122_Get_bed.sh  
GS20240122_KDii_process_input_files.sh  
GS20240122_KDii_extract_reads.sh  
GS20240122_KDii_extract_reads_ii.sh  
GS20240122_Get_bamlist.sh  
GS20240122_KDii_run_GS.sh  
GS20240122_MAP_CIs_from_topbottom_KD_ii.sh  

Scripts to visualize results  
GS20240122_KDii_visualization.R  

Scripts for differential expression  
GS20240122_featureCounts.sh  
GS20240122_DEseq2.Rmd  

Retained introns GRAND-SLAM control runs  
GS20240122_Get_retained_intron_gtf.sh  
GS20240122_Get_RI_gedi_oml.sh  
  
scr pc parameter sweep control runs  
GS20240305_Get_pc_files_folders.sh  
GS20240305_KDii_run_GS.sh  
Comparison_GS_KDii_20240305.ipynb  

Final run for MS  
GS20240309_Get_pc_files_folders.sh  
GS20240309_KDii_run_GS.sh  


# 10_lasso
First round of feature selection using individual feature classes with lasso regression  
ML_20240122_subcell_feat_select1.sh  
which calls: ML_20240122_subcell_feat_select1.py

Second round feature selection using the union of features selected in round 1  
ML_20240122_subcell_feat_select2.sh  
which calls: ML_20240122_subcell_feat_select2.py

Optimal hyperparameter identification, visualization and final file output  
ML_20240122_subcell_analysis.ipynb  
