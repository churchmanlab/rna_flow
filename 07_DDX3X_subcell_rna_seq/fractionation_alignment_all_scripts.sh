# Churchman fractionation alignment

# They use refseq GTF (annotations) https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_000001405.40/ and the ENSEMBL GRCh38 https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/ 
# They also uploaded this GTF to basespace which is from the link so I am using their uploaded version.
# Homo_sapiens.GRCh38.104.gtf.zip
# 47.4 MB Download

# according to STAR manual https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf ENSEMBL Fasta and UCSC Fastas have diff chromosome naming convertions (UCSCS uses chr1,2, ensembl just uses 1,2,) so I have to re-download the ensembl Fasta of GRCh38

wget https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

gzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa



# With --quantMode GeneCounts option STAR will count number reads per gene while mapping. A
# read is counted if it overlaps (1nt or more) one and only one gene. Both ends of the paired-end
# read are checked for overlaps. The counts coincide with those produced by htseq-count with default
# parameters. This option requires annotations in GTF format (i.e. gene id tag for each exon) specified
# in --sjdbGTFfile at the genome generation step or at the mapping step provided in option. STAR
# outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to
# different strandedness options:
# column 1: gene ID
# column 2: counts for unstranded RNA-seq
# column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)


# 2 x 151bp PE Novaseq S4 reads so set sjdbOverhang to 151-1 = 150

###################################
### generate_STAR_indices_150.sh ###
###################################
#!/bin/env bash
module load CBI star/2.7.5a

set -x
STAR --runThreadN $NSLOTS \
--runMode genomeGenerate \
--genomeDir /wynton/home/floorlab/axu/references/STAR_indices_150_ENSEMBL-Refseq/ \
--genomeFastaFiles /wynton/home/floorlab/axu/references/ENSEMBL_Refseq_references/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile /wynton/home/floorlab/axu/references/ENSEMBL_Refseq_references/Homo_sapiens.GRCh38.104.gtf \
--sjdbOverhang 150 

qstat -j $JOB_ID
# end

qsub -cwd -pe smp 8 -l mem_free=8G -l h_rt=10:00:00 generate_STAR_indices_150.sh
# 3 hrs 6 core 8G ram





######################
### makefolders.sh ###
######################
sample_name=$1

echo "$sample_name"
mkdir "$sample_name"
mv "$sample_name"_R1.fastq.gz "$sample_name"
mv "$sample_name"_R2.fastq.gz "$sample_name"
# end script

cat samplename.txt | xargs -n 1 sh -c 'bash makefolders.sh $0'


########################
### RNA_alignment.sh ###
########################

# Uses STAR to make alignments
#!/bin/env bash
set -x

module load CBI star/2.7.5a
module load CBI samtools/1.10
# module load CBI cutadapt/4.1

sample_name=$1
cd "$sample_name"

# Adaptor trimming sequences:
# The NEBNext libraries for Illumina resemble TruSeq libraries and can be trimmed similar to TruSeq:
# Adaptor Read 1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# Adaptor Read 2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
# cutadapt \
#     -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
#     -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     --cores $NSLOTS \
#     -o "$sample_name"_trimmed_R1.fastq.gz -p "$sample_name"_trimmed_R2.fastq.gz \
#     "$sample_name"_R1.fastq.gz "$sample_name"_R2.fastq.gz
# weird STAR error after using CutAdapt. might be trimming reads to 0 or something. will switch to using STAR's built in method
# EXITING because of FATAL ERROR in reads input: short read sequence line: 0
# Read Name=@A01102:540:HK2WNDSX5:2:1101:2474:3521
# Read Sequence====
# DEF_readNameLengthMax=50000
# DEF_readSeqLengthMax=650

# lorenzo's paramters
STAR --genomeDir /wynton/home/floorlab/axu/references/STAR_indices_150_ENSEMBL-Refseq/ \
--readFilesCommand zcat \
--readFilesIn "$sample_name"_R1.fastq.gz "$sample_name"_R2.fastq.gz \
--runThreadN $NSLOTS \
--outFilterMismatchNmax 3 \
--outFilterMultimapNmax 50 \
--chimScoreSeparation 10 \
--chimScoreMin 20 \
--chimSegmentMin 15 \
--outSAMattributes NH HI AS nM NM MD XS \
--outFilterIntronMotifs RemoveNoncanonicalUnannotated \
--alignSJoverhangMin 500 \
--outFileNamePrefix $sample_name"_" \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--outSAMmultNmax 1 \
--outMultimapperOrder Random \
--limitBAMsortRAM 19000000000 \
--alignEndsType Local \
--quantMode GeneCounts \
--clip3pAdapterSeq "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" \
--clip3pAdapterMMp 0.1 0.1 

samtools index -@ $NSLOTS "$sample_name"_Aligned.sortedByCoord.out.bam

qstat -j $JOB_ID
# end script

# submit parallel jobs to wynton
cat samplename.txt | xargs -n 1 sh -c 'qsub -cwd -pe smp 4 -l mem_free=4G -l h_rt=6:00:00 RNA_alignment.sh $0'


# organizing
mkdir STAR_logs
find . -type f -name "*.final.out" -exec cp {} STAR_logs \;

mkdir ReadsPerGene_All
find . -type f -name "*ReadsPerGene.out.tab" -exec cp {} ReadsPerGene_All \;

mkdir all_bams
find . -type f -name "*.bam*" -exec cp {} all_bams \;

scp -r axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu-dup/AXSNF10_cutadapted/STAR_logs/ .
scp -r axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu-dup/AXSNF10_cutadapted/ReadsPerGene_All/ .

scp axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu/AXSNF10/all_bams/flag_1_H9_Aligned.sortedByCoord.out.bam.bai .
scp axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu/AXSNF10/all_bams/flag_1_H9_Aligned.sortedByCoord.out.bam .


scp axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu/AXSNF10/all_bams/K562_DDX3-R376C_Nuc_rep1_Aligned.sortedByCoord.out.bam.bai .
scp axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu/AXSNF10/all_bams/K562_DDX3-R376C_Nuc_rep1_Aligned.sortedByCoord.out.bam .

scp axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu/AXSNF10/all_bams/K562_DDX3-R376C_Cyto_rep1_Aligned.sortedByCoord.out.bam.bai .
scp axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu/AXSNF10/all_bams/K562_DDX3-R376C_Cyto_rep1_Aligned.sortedByCoord.out.bam .




# run multiqc on STAR logs
cd /Users/albert/Dropbox/Floor_Lab/Analysis/AXSNF10_Churchman_collab/STAR_logs/
/Users/albert/Library/Python/3.8/bin/multiqc .


