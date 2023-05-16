pip3 install --upgrade autopep8 --user
albert@Alberts-MacBook-Pro-3 ~ % cd /Users/albert/Library/Python/3.8/lib/python/site-packages
albert@Alberts-MacBook-Pro-3 site-packages % chmod +x ./autopep8.py

# have to reformat fastqgz_to_counts.py using autopep8 to get it to run on mac os python3. some tab/space indention issue. This formatting different than the one from Max's github.
albert@Alberts-MacBook-Pro-3 ScreenProcessing % python3 /Users/albert/Library/Python/3.8/lib/python/site-packages/autopep8.py -i /Users/albert/Dropbox/Floor_Lab/Analysis/Programs/ScreenProcessing/fastqgz_to_counts.py

#fastqgz to counts for sgDDX3 (using reformatted fastqgz_to_counts.py)

python3 /Users/albert/Dropbox/Floor_Lab/Analysis/Programs/ScreenProcessing/fastqgz_to_counts.py -p 2 --trim_start 1 --trim_end 29 /Users/albert/Dropbox/Gilbert\ lab\ rotation/ScreenProcessing-master/library_reference/CRISPRi_v2_human.trim_1_29_forward.fa /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/output_python3 /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/*.fastq.gz

# Library file loaded successfully:
# 	2.09E+05 elements (2.06E+05 unique sequences)	28bp reads expected
# Processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/sgDDX3_T0_rep1.fastq.gz
# Processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/sgDDX3_T0_rep2.fastq.gz
# Done processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/sgDDX3_T0_rep1.fastq.gz
# Processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/sgDDX3_Tf_rep1.fastq.gz
# Done processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/sgDDX3_T0_rep2.fastq.gz
# Processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/sgDDX3_Tf_rep2.fastq.gz
# Done processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/sgDDX3_Tf_rep1.fastq.gz
# Done processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/sgDDX3_Tf_rep2.fastq.gz
# /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/output_python3/count_files/sgDDX3_T0_rep1_CRISPRi_v2_human.trim_1_29_forward.fa.counts:
# 	6.89E+07 reads	5.94E+07 aligning (86.23%)
# /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/output_python3/count_files/sgDDX3_T0_rep2_CRISPRi_v2_human.trim_1_29_forward.fa.counts:
# 	8.64E+07 reads	7.40E+07 aligning (85.66%)
# /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/output_python3/count_files/sgDDX3_Tf_rep1_CRISPRi_v2_human.trim_1_29_forward.fa.counts:
# 	7.25E+07 reads	6.14E+07 aligning (84.62%)
# /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgDDX3/output_python3/count_files/sgDDX3_Tf_rep2_CRISPRi_v2_human.trim_1_29_forward.fa.counts:
# 	7.68E+07 reads	6.55E+07 aligning (85.31%)
# Done processing all sequencing files



python3 /Users/albert/Dropbox/Floor_Lab/Analysis/Programs/ScreenProcessing/fastqgz_to_counts.py -p 2 --trim_start 1 --trim_end 29 /Users/albert/Dropbox/Gilbert\ lab\ rotation/ScreenProcessing-master/library_reference/CRISPRi_v2_human.trim_1_29_forward.fa /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/output_python3 /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/*.fastq.gz

# Library file loaded successfully:
# 	2.09E+05 elements (2.06E+05 unique sequences)	28bp reads expected
# Processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/sgGal4_T0_rep1.fastq.gz
# Processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/sgGal4_T0_rep2.fastq.gz
# Done processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/sgGal4_T0_rep1.fastq.gz
# Processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/sgGal4_Tf_rep1.fastq.gz
# Done processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/sgGal4_T0_rep2.fastq.gz
# Processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/sgGal4_Tf_rep2.fastq.gz
# Done processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/sgGal4_Tf_rep1.fastq.gz
# Done processing /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/sgGal4_Tf_rep2.fastq.gz
# /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/output_python3/count_files/sgGal4_T0_rep1_CRISPRi_v2_human.trim_1_29_forward.fa.counts:
# 	4.89E+07 reads	4.18E+07 aligning (85.48%)
# /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/output_python3/count_files/sgGal4_T0_rep2_CRISPRi_v2_human.trim_1_29_forward.fa.counts:
# 	6.56E+07 reads	5.61E+07 aligning (85.54%)
# /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/output_python3/count_files/sgGal4_Tf_rep1_CRISPRi_v2_human.trim_1_29_forward.fa.counts:
# 	7.59E+07 reads	6.54E+07 aligning (86.11%)
# /Users/albert/Desktop/Sequencing/AXSNF03-P/210617_K00153_0911_BHLGNJBBXY_SE50/AXSNF03-P/concatenated_reads/sgGal4/output_python3/count_files/sgGal4_Tf_rep2_CRISPRi_v2_human.trim_1_29_forward.fa.counts:
# 	6.48E+07 reads	5.57E+07 aligning (85.98%)

# had to install configparser and pandas - in the future maybe switch to running everything in python 3
pip3 install configparser
pip3 install pandas
pip3 install scipy
pip3 install matplotlib

# had to change applyMW function in process_experiments.py (from Marco Jost) to accomodate new SciPy version 1.7.1 See process_experiments_new.py for old and new version of function


# screen processing sgDDX3
python3 /Users/albert/Dropbox/Floor_Lab/Analysis/Programs/ScreenProcessing/process_experiments_new.py /Users/albert/Dropbox/Floor_Lab/Analysis/AXSNF03/AXSNF03_sgDDX3_config_file_python3.txt /Users/albert/Dropbox/Gilbert\ lab\ rotation/ScreenProcessing-master/library_tables
# No growth values--all phenotypes will be reported as log2enrichments

# Accessing library information
# Loading counts data
# Merging experiment counts split across lanes/indexes
# -generating sgRNA read count histograms
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/000_fig_counts_hist.png
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/001_fig_counts_hist.png
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/002_fig_counts_hist.png
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/003_fig_counts_hist.png
# Computing sgRNA phenotype scores
# -generating phenotype histograms and scatter plots
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/004_fig_counts_scatter.png
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/005_fig_phenotype_hist.png
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/006_fig_sgRNAs_passing_filter_hist.png
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/007_fig_counts_scatter.png
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/008_fig_phenotype_hist.png
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/009_fig_sgRNAs_passing_filter_hist.png
# Averaging replicates
# -generating replicate phenotype histograms and scatter plots
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/010_fig_phenotype_scatter.png
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/011_fig_phenotype_hist.png
# Generating a pseudogene distribution from negative controls
# Computing gene scores
# --calculate_ave
# --calculate_mw
# Collapsing transcript scores to gene scores
# AXSNF03_sgDDX3_results/AXSNF03_sgDDX3_results_plots/012_fig_volcano_plot.png
# Done!


# screen processing sgGal4
python3 /Users/albert/Dropbox/Floor_Lab/Analysis/Programs/ScreenProcessing/process_experiments_new.py /Users/albert/Dropbox/Floor_Lab/Analysis/AXSNF03/AXSNF03_sgGal4_config_file_python3.txt /Users/albert/Dropbox/Gilbert\ lab\ rotation/ScreenProcessing-master/library_tables


# No growth values--all phenotypes will be reported as log2enrichments

# Accessing library information
# Loading counts data
# Merging experiment counts split across lanes/indexes
# -generating sgRNA read count histograms
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/000_fig_counts_hist.png
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/001_fig_counts_hist.png
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/002_fig_counts_hist.png
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/003_fig_counts_hist.png
# Computing sgRNA phenotype scores
# -generating phenotype histograms and scatter plots
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/004_fig_counts_scatter.png
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/005_fig_phenotype_hist.png
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/006_fig_sgRNAs_passing_filter_hist.png
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/007_fig_counts_scatter.png
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/008_fig_phenotype_hist.png
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/009_fig_sgRNAs_passing_filter_hist.png
# Averaging replicates
# -generating replicate phenotype histograms and scatter plots
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/010_fig_phenotype_scatter.png
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/011_fig_phenotype_hist.png
# Generating a pseudogene distribution from negative controls
# Computing gene scores
# --calculate_ave
# --calculate_mw
# Collapsing transcript scores to gene scores
# AXSNF03_sgGal4_results/AXSNF03_sgGal4_results_plots/012_fig_volcano_plot.png
# Done!




