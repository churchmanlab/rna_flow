#!/bin/bash
#Robert Ietswaart 
#Modify the existing SNP-masked gtf to genes only.


date=20240126
baseDir=/n/groups/churchman/ri23/bseq/GS${date}
  
module load gcc/9.2.0
module load java/jdk-1.8u112
module load bedops/2.4.30
module load bedtools/2.30.0
module load samtools/1.15.1

cd ${baseDir}
mkdir -p ${baseDir}/LogErr

input_gtf=K562_ensGRCh38_MTmod_dm6_ercc_cat.gtf
output_gtf=K562_ensGRCh38_MTmod_dm6_ercc_cat_genes_only.gtf

# #0 Copy starting gtf with all feature to local folder
cp /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/$input_gtf ./$input_gtf

#1 output only the gene features from a GTF file and label gene as exon so GS 
#GRAND-SLAM gtf requirements: https://github.com/erhard-lab/gedi/wiki/Preparing-genomes
#Note that GS require to have CDS and start/stop codons in the file for protein coding genes

awk 'BEGIN{FS=OFS="\t"} $1 ~ /^#!/ || $3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon" {print; next} $3 == "gene" {match($9, /gene_id "([^"]+)"/, arr); $3="exon"; $9="gene_id \""arr[1]"\"; transcript_id \"t_"arr[1]"\";"; print}' $input_gtf > $output_gtf



#OLD:
#Gives error: transcripts.fasta is empty (due to missing transcript_id)
# awk 'BEGIN{FS=OFS="\t"} $1 ~ /^#!/ || $3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon" {print; next} $3 == "gene" {$3="exon"; print}' $input_gtf > $output_gtf

# #Gives error: transcript_id not consecutive in file
# awk 'BEGIN{FS=OFS="\t"} $1 ~ /^#!/ || $3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon" {print; next} $3 == "gene" {$3="exon"; $9=$9" transcript_id \"t_"$1"_"$4"_"$5"\";"; print}' $input_gtf > $output_gtf

