#!/bin/bash

#SBATCH -c 3
#SBATCH -t 0-24:00
#SBATCH -p medium
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mtcouvi@gmail.com
#SBATCH -e logs/CallSNPs_ModifyGenome_%j.err
#SBATCH -o logs/CallSNPs_ModifyGenome_%j.log


### Written 2020 by M. Couvillion
### This script will call snps with bcftools using one or more bam files, separate the resulting .vcf file into INDELs and non-INDELs, and modify the reference fasta based on the non-INDEL snvs. A custom python script is then used to further split the vcf into loci that have > 75% of coverage as a specific alternate base (these will be changed in the reference to the alternate base), and those that have < 75 % an alternate base (these will be changed to N in the reference). This first reference modification is done with gatk FastaAlternateReferenceMaker. Next a pass filter column is added to the INDEL vcf with another custom python script and another round of reference modification is done from that using rf2m which modifies the matching gtf annotation file in parallel. Finally, matching fasta index, dictionary, and chrom.sizes files are made.
### Updated 8/23 to use only primary alignments, also increased quality for SNP and INDEL selection ** new ParseAndFilterVCF.py and ParseAndPASSlabelVCF_v2.py

### USE
# Modify this script: Naming parameters, the list of bams, the filtering for primary alignments, and the list of bams in the bcftools mpileup command, then:
# sbatch ../Scripts/CallSNPs_ModifyGenome.sh 
# run from directory with .bam files





# Load required module
module load picard/2.8.0
module load samtools/1.9
module load bcftools/1.9
module load gatk/4.0.0.0
# module load snpEff/4.3g

# helpful tips: https://www.biostars.org/p/327558/

# Set naming paramters
Celline="K562" # Hela K562 mouseNIH3T3 HEK293T
MapMethod="t5MMinformed5" #t5MMinformed5 MMinformed5
# Set quality cutoff
qual1=20
qual2=30
AFcutoff=0.75

# Set path to sequence directory
path="/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/"

# Set paths to fasta
infasta='/n/groups/churchman/mc348/TimelapseSeq/SeqFiles/ensGRCh38_h_MT_ncRNAs_allERCC_merge.formatted.fasta' 


out1fasta=${path}ensGRCh38_h_MT_ncRNAs_allERCC_merge.${Celline}snpCor.fasta
# mm10_dm6_ercc_cat.snpCor.fasta
# ensGRCh38_dm6_ercc_cat.snpCor.fasta 
# ensGRCh38_h_MT_ncRNAs_allERCC_merge.snpCor.fasta
out2fasta=${path}${Celline}_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta
# ${Celline}_mm10_dm6_ercc_cat.fasta
# ${Celline}_ensGRCh38_dm6_ercc_cat.fasta
# ${Celline}_ensGRCh38_h_MT_ncRNAs_allERCC_merge.fasta

# Set paths to gtf
ingtf=${path}ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf
# ingtf=${path}Hela_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf  # Use this for new human genomes since it has anti genes, 7S, etc 08162023: No, the original version of this is not quite right
# mm10_MTmod_dm6_ercc_cat.gtf 
# ensGRCh38_fusedMTgenes_dm6_ercc_cat.gtf
# ensGRCh38_h_MT_ncRNAs_allERCC_merge_fusedMTgenes2.gtf
outgtf=${path}${Celline}_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf
# ${Celline}_mm10_MTmod_dm6_ercc_cat.gtf 
# ${Celline}_ensGRCh38_MTmod_dm6_ercc_cat.gtf  
# ${Celline}_ensGRCh38_h_MT_ncRNAs_allERCC_merge_MTmod.gtf

# Need to make sure there is a fasta reference and dictionary for this file in the same directory, e.g.
# faidx reference:
rm ${infasta}.fai
samtools faidx $infasta
# # fasta dictionary
rm ${infasta/fasta/dict}
java -jar $PICARD/picard-2.8.0.jar CreateSequenceDictionary R=$infasta O=${infasta/fasta/dict}

# Set path to input files  ** Make sure these have been mapped to infasta above **

# For Hela, repeat 06162023
# inbam1='/n/groups/churchman/mc348/TimelapseSeq/Hela_forGenome/TL1_0m_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
# inbam2='/n/groups/churchman/mc348/TimelapseSeq/Hela_forGenome/TL3_0m_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
# inbam3='/n/groups/churchman/mc348/TimelapseSeq/Hela_forGenome/TL4_0m_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
# inbam4='/n/groups/churchman/mc348/TimelapseSeq/Hela_forGenome/TL5_0m_tot_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'




# For HEK
# inbam1='/n/groups/churchman/mc348/TimelapseSeq/HEK293_forGenome/TL10_0m_WT_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
# inbam2='/n/groups/churchman/mc348/TimelapseSeq/HEK293_forGenome/TL12_0m_WT_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
# inbam3='/n/groups/churchman/mc348/TimelapseSeq/HEK293_forGenome/TL14_0m_WT_A_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
# inbam4='/n/groups/churchman/mc348/TimelapseSeq/HEK293_forGenome/TL14_0m_WT_B_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'


# For K562
inbam1='/n/groups/churchman/mc348/TimelapseSeq/K562_forGenome/T1-tot_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
inbam2='/n/groups/churchman/mc348/TimelapseSeq/K562_forGenome/T1-poly_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
inbam3='/n/groups/churchman/mc348/TimelapseSeq/K562_forGenome/T1-chr_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
inbam4='/n/groups/churchman/mc348/TimelapseSeq/K562_forGenome/U1-tot_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
inbam5='/n/groups/churchman/mc348/TimelapseSeq/K562_forGenome/U1-poly_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'
inbam6='/n/groups/churchman/mc348/TimelapseSeq/K562_forGenome/U1-chr_t5MMinformed5_Aligned.sortedByCoord.noSpike.bam'

# Remove non-primary alignments and singletons
echo 'Filtering out non-primary alignments'
echo 'Working on bam 1'
inbam1p=${inbam1/bam/primary.bam}
samtools view -@ 3 -F 0x8 -F 0x100 -o $inbam1p $inbam1
echo 'Working on bam 2'
inbam2p=${inbam2/bam/primary.bam}
samtools view -@ 3 -F 0x8 -F 0x100 -o $inbam2p $inbam2
echo 'Working on bam 3'
inbam3p=${inbam3/bam/primary.bam}
samtools view -@ 3 -F 0x8 -F 0x100 -o $inbam3p $inbam3
echo 'Working on bam 4'
inbam4p=${inbam4/bam/primary.bam}
samtools view -@ 3 -F 0x8 -F 0x100 -o $inbam4p $inbam4
echo 'Working on bam 5'
inbam5p=${inbam5/bam/primary.bam}
samtools view -@ 3 -F 0x8 -F 0x100 -o $inbam5p $inbam5
echo 'Working on bam 6'
inbam6p=${inbam6/bam/primary.bam}
samtools view -@ 3 -F 0x8 -F 0x100 -o $inbam6p $inbam6
echo 'Finished filtering bams'


# Call snps 
# This only calls loci with >1 read
# mpileup generates genotype likelihoods at each genomic position with coverage
echo 'Starting variant calling with bcftools'
bcftools mpileup -Ou -f $infasta $inbam1p $inbam2p $inbam3p $inbam4p $inbam5p $inbam6p | bcftools call -mv -Ov -o ${Celline}_${MapMethod}_snpcalls.vcf # Did not finish in 12 hours with all 8 bams
echo 'Finished variant calling'
# # 
# echo Finished making snpcalls.vcf
# 
# # make separate files for INDELs and non-INDELS
grep -v 'INDEL' ${Celline}_${MapMethod}_snpcalls.vcf > ${Celline}_${MapMethod}_snpcalls.noINDEL.vcf
grep -e '^#' -e 'INDEL' ${Celline}_${MapMethod}_snpcalls.vcf > ${Celline}_${MapMethod}_snpcalls.INDELonly.vcf
# # make separate files for those that get the alternate allele and those that get N (must have at least 5 reads to get alternate allele)
python /n/groups/churchman/mc348/TimelapseSeq/Scripts/ParseAndFilterVCF.py ${Celline}_${MapMethod}_snpcalls.noINDEL.vcf $AFcutoff $qual1
# # # # # 
# # # # # index noINDEL vcfs for next step
/n/groups/churchman/mc348/programs/gatk-4.1.9.0/gatk IndexFeatureFile -I ${Celline}_${MapMethod}_snpcalls.noINDEL.AFatleast${AFcutoff}.vcf
/n/groups/churchman/mc348/programs/gatk-4.1.9.0/gatk IndexFeatureFile -I ${Celline}_${MapMethod}_snpcalls.noINDEL.AFbelow${AFcutoff}.vcf
# # # # 
# # # # # First fasta modification without considering INDELs (takes ~30min) 
# # # # This uses gatk to modify fasta with Ns for SNPs with ambiguous alternate allele and the alternate base for SNPs where it's unambiguous
/n/groups/churchman/mc348/programs/gatk-4.1.9.0/gatk FastaAlternateReferenceMaker -R $infasta -O $out1fasta -V ${Celline}_${MapMethod}_snpcalls.noINDEL.AFatleast${AFcutoff}.vcf --snp-mask ${Celline}_${MapMethod}_snpcalls.noINDEL.AFbelow${AFcutoff}.vcf --snp-mask-priority true 
# # # 
# # # # Have to delete the number added in front of each fasta header, e.g. 
# # # # >1 1:1-249250621
# # # # >2 10:1-135534747
# # # # >3 h_MT:1-16023
sed -i 's/>.* />/' $out1fasta
sed -i 's/:/ /' $out1fasta
# 
# # # Remove non-matched index and dictionary files 
rm ${out1fasta/fasta/dict}
rm ${out1fasta/fasta/fasta.fai}

# Next use rf2m to deal with INDELs - this makes a new fasta with matched gtf since coordinates will change when INDELs are allowed

# First filter by indel quality, proximity to each other, and length, and add PASS to FILTER column for rf2m to read line (number is quality filter cutoff for PASS, output filename ends in '.QUALfilt${qual}.vcf')
python /n/groups/churchman/mc348/TimelapseSeq/Scripts/ParseAndPASSlabelVCF_v2.py ${Celline}_${MapMethod}_snpcalls.INDELonly.vcf $qual2 ${AFcutoff}



# Using rf2m to make a new fasta and gtf together, ***including indels***
module load perl/5.30.0

/n/groups/churchman/mc348/programs/rf2m-master/format_vcf.pl --vcf ${Celline}_${MapMethod}_snpcalls.INDELonly.QUALfilt${qual2}.vcf

/n/groups/churchman/mc348/programs/rf2m-master/genome_creator_corrected.pl --genome $out1fasta --outGenome $out2fasta      ### Had to modify the genome_creator.pl script to 
# 1. read gencode fasta (originally written for ensembl fasta:  $header =~ m/^(.+) dna/; >  $header =~ m/^(.+) 1-/;
# 2. Exit loop if a base has already been changed in the gatk step above, otherwise this creates an offset in the gtf later (line 71: next > last)


/n/groups/churchman/mc348/programs/rf2m-master/gtf_creator.pl --gtf $ingtf --outGTF $outgtf



# Make new fasta index and dictionary
# Make fasta index
samtools faidx $out2fasta
# Make sequence dictionary
rm ${out2fasta/fasta/dict}
java -jar $PICARD/picard-2.8.0.jar CreateSequenceDictionary R=$out2fasta O=${out2fasta/fasta/dict}

# Make new chrom.sizes for modified fasta:
cut -f1,2 ${out2fasta}.fai > ${path}${Celline}.chrom.sizes


# to make matching ENS to GeneName file:
if [ "${Celline}" = "NIH3T3" ]
then
script="/n/groups/churchman/mc348/TimelapseSeq/Scripts/GENCODE_gtf2ENSMUSGtoGeneNameAndExonLengths.py"
else
script="/n/groups/churchman/mc348/TimelapseSeq/Scripts/GENCODE_gtf2ENStoGeneNameAndExonLengths.py"
fi


output1=${path}${Celline}_ENStoGeneName.txt
output2=${path}${Celline}_ExonLengths.txt

rm $output1
rm $output2
python ${script} -i $outgtf -o $output1 -O $output2

# Make .bed file
if [ "${Celline}" = "NIH3T3" ]
then
species="mouse"
else																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																
species="human"
fi


script="/n/groups/churchman/mc348/TimelapseSeq/Scripts/GENCODE_gtf2bed.py"
bedname=${outgtf/gtf/bed}
rm $bedname
python $script -i $outgtf -o $bedname -g $species

# Make .bed file with outside-most coordinates
script="/n/groups/churchman/mc348/TimelapseSeq/Scripts/ChooseMaxTxpt.R"
rm ${bedname/.bed/_outsideMostCoords.bed}																																																																				
Rscript $script $bedname ${bedname/.bed/_outsideMostCoords.bed}

# To make fasta and gtf for loading on IGV:
fastagz=${out2fasta/fasta/fasta.gz}
gzip -c $out2fasta > $fastagz
sbatch /n/groups/churchman/mc348/TimelapseSeq/Scripts/index-gtf.sh $outgtf


# Clean up
# rm *primary.bam
rm $out1fasta
# rm *vcf
