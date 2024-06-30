#!/bin/bash
#Robert Ietswaart 
#Modify the existing SNP-masked gtf to filter for (all/no exons) from transcript isoforms containing retained introns.
#Source: GS20240122_Get_gene_level_gtf.sh


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
output_gtf1=K562_ensGRCh38_MTmod_dm6_ercc_cat_protein_coding.gtf
output_gtf2=K562_ensGRCh38_MTmod_dm6_ercc_cat_no_retained_introns.gtf
output_gtf3=K562_ensGRCh38_MTmod_dm6_ercc_cat_retained_introns.gtf
output_gtf4=K562_ensGRCh38_MTmod_dm6_ercc_cat_retained_introns_with_exons.gtf


# #0 Copy starting gtf with all feature to local folder: already done in GS20240122_Get_gene_level_gtf.sh 
# cp /n/groups/churchman/mc348/TimelapseSeq/SeqFiles/$input_gtf ./$input_gtf

#1 filter for exons originating from PC isoforms only.
#Note that GS require to have CDS and start/stop codons in the file for protein coding genes
awk 'BEGIN{FS=OFS="\t"} $1 ~ /^#!/ || $3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon" || $3 == "exon" && /transcript_biotype "protein_coding"/' $input_gtf > $output_gtf1

#2 "minimally invasive" procedure that only filters out exons corresponding to RI isoforms. Shared exons will still be present as part of PC (and other type) isoforms
awk '!/transcript_biotype "retained_intron"/' $input_gtf > $output_gtf2

#3 filter for "exons" which are exclusively present in the RI isoform, so these "exons" should represent the RIs 
awk 'BEGIN{FS=OFS="\t"} $1 ~ /^#!/ || $3 == "CDS" || $3 == "start_codon" || $3 == "stop_codon" || $3 == "exon" && /transcript_biotype "retained_intron"/' $input_gtf > retained_introns_temp.gtf
awk 'BEGIN{FS=OFS="\t"} $3 == "exon" && !/transcript_biotype "retained_intron"/' $input_gtf > not_retained_introns_temp.gtf
cut -f 1,3,4,5,7 retained_introns_temp.gtf > retained_introns.coords
cut -f 1,3,4,5,7 not_retained_introns_temp.gtf > not_retained_introns.coords
grep -F -f not_retained_introns.coords retained_introns.coords > shared_exons.coords
awk 'BEGIN{FS="\t"} NR==FNR {a[$0]; next} {c=$1FS$3FS$4FS$5FS$7; if (!(c in a)) print }' shared_exons.coords retained_introns_temp.gtf > $output_gtf3

# #4 filter for "exons" which are both non- and exclusively present in the RI isoform
cp retained_introns_temp.gtf $output_gtf4

#rm not_retained_introns_temp.gtf
#rm retained_introns_temp.gtf