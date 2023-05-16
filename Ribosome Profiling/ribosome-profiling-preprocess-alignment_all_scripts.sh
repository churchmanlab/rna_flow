cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgDDX3_RP_rep1_d0_S3_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgDDX3_RP_rep1_d0_S72_L006_R1_001.fastq.gz > sgDDX3_RP_rep1_d0.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgDDX3_RP_rep1_d14_S11_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgDDX3_RP_rep1_d14_S80_L006_R1_001.fastq.gz > sgDDX3_RP_rep1_d14.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgDDX3_RP_rep1_d3_S7_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgDDX3_RP_rep1_d3_S76_L006_R1_001.fastq.gz > sgDDX3_RP_rep1_d3.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgDDX3_RP_rep2_d0_S4_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgDDX3_RP_rep2_d0_S73_L006_R1_001.fastq.gz > sgDDX3_RP_rep2_d0.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgDDX3_RP_rep2_d14_S12_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgDDX3_RP_rep2_d14_S81_L006_R1_001.fastq.gz > sgDDX3_RP_rep2_d14.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgDDX3_RP_rep2_d3_S8_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgDDX3_RP_rep2_d3_S77_L006_R1_001.fastq.gz > sgDDX3_RP_rep2_d3.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgGal4_RP_rep1_d0_S1_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgGal4_RP_rep1_d0_S70_L006_R1_001.fastq.gz > sgGal4_RP_rep1_d0.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgGal4_RP_rep1_d14_S9_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgGal4_RP_rep1_d14_S78_L006_R1_001.fastq.gz > sgGal4_RP_rep1_d14.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgGal4_RP_rep1_d3_S5_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgGal4_RP_rep1_d3_S74_L006_R1_001.fastq.gz > sgGal4_RP_rep1_d3.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgGal4_RP_rep2_d0_S2_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgGal4_RP_rep2_d0_S71_L006_R1_001.fastq.gz > sgGal4_RP_rep2_d0.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgGal4_RP_rep2_d14_S10_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgGal4_RP_rep2_d14_S79_L006_R1_001.fastq.gz > sgGal4_RP_rep2_d14.fastq.gz
cat /bfd/axu/raw_sequencing/AXSNF07_210727/sgGal4_RP_rep2_d3_S6_L001_R1_001.fastq.gz /bfd/axu/raw_sequencing/AXSNF07_210804/sgGal4_RP_rep2_d3_S75_L006_R1_001.fastq.gz > sgGal4_RP_rep2_d3.fastq.gz

# download all to wynton
scp -r axu@brubeck.ucsf.edu:/bfd/axu/raw_sequencing/AXSNF07_concat/ .

# extract all files in a folder with .gz
gzip -d *.gz


###################
### gzip_all.sh ###
###################
set -x

ls *.gz
gzip -d *.gz # extract all files in a folder with .gz

qstat -j $JOB_ID
## end script ##

#submit job
qsub -cwd -pe smp 4 -l mem_free=4G -l h_rt=3:00:00 /wynton/scratch/axu/gzip_all.sh # ~20 min with 4 cores, 4G



# make folder for every file with name (before .) and move file into folder
for file in *; do
  if [[ -f "$file" ]]; then
    mkdir "${file%.*}"
    mv "$file" "${file%.*}"
  fi
done


# install miniconda3 for cutadapt
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

/wynton/home/floorlab/axu/miniconda3/bin/conda config --add channels defaults
/wynton/home/floorlab/axu/miniconda3/bin/conda config --add channels bioconda
/wynton/home/floorlab/axu/miniconda3/bin/conda config --add channels conda-forge
/wynton/home/floorlab/axu/miniconda3/bin/conda install -c bioconda cutadapt


# install fastx toolkit
wget http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2
tar xjvf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2


# make a samplenames txt
cat samplenames.txt
sgDDX3_RP_rep1_d0
sgDDX3_RP_rep1_d14
sgDDX3_RP_rep1_d3
sgDDX3_RP_rep2_d0
sgDDX3_RP_rep2_d14
sgDDX3_RP_rep2_d3
sgGal4_RP_rep1_d0
sgGal4_RP_rep1_d14
sgGal4_RP_rep1_d3
sgGal4_RP_rep2_d0
sgGal4_RP_rep2_d14
sgGal4_RP_rep2_d3


# added these lines to contaminant_rnaseq_new.fa to make contaminant_rnaseq_new_sgRNAs.fa
>sgGal4_pAX71
GAACGACTAGTTAGGCGTGTAGTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTC
>sgDDX3_pAX90
GACCGCGAAGGCCCTCTCACGTTTAAGAGCTAAGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTC




#######################
### bowtie_build.sh ###
#######################
# builds bowtie index for contaminant RNA file

#!/bin/evn bash
module load CBI bowtie2/2.4.1 
bowtie2-build --threads $NSLOTS contaminant_rnaseq_new_sgRNAs.fa contaminant_rnaseq_new_sgRNAs

qstat -j $JOB_ID
## end script ##




############################################
### RP_pre_alignment_withsgRNA_Filter.sh ###
############################################

module load CBI bowtie2/2.4.1
set -x

sample_name=$1
fastq="${sample_name}.fastq"
cd $sample_name

# remove sgDDX3 sequence if applicable
# literal= is where I see on IGV the most reads - corresponds to most of the protospacer
/wynton/home/floorlab/axu/bbmap/bbduk.sh \
in=$fastq \
out=${sample_name}_sgRNAfilter.fastq \
outm=sgRNA_contam.fastq \
literal=CGCGAAGGCCCTCTCA \
k=14 \
hdist=0 \
stats=${sample_name}_bbduk_stats.txt

# convert from fastq to fasta
/wynton/home/floorlab/axu/programs/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/bin/fastq_to_fasta -Q33 -i ${sample_name}_sgRNAfilter.fastq -o ribo.fasta

# remove adapter sequence (?)
/wynton/home/floorlab/axu/miniconda3/bin/cutadapt -m 18 --discard-untrimmed -a TAGACAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o trim_cut.fasta ribo.fasta > cutadapt.log

# collapse reads by UMI
/wynton/home/floorlab/axu/programs/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/bin/fastx_collapser -i trim_cut.fasta -o coll_trim_cut.fasta

# remove UMI
# filename out: *-noUMI.fasta 
# coll_trim_cut-noUMI.fasta
perl /wynton/home/floorlab/axu/scripts/remove5N_3Nfasta.pl coll_trim_cut.fasta  1> rem44_indexes.log

# align against repeat RNAs
bowtie2 -f -k 1 --local -p 2 -x /wynton/home/floorlab/axu/references/contaminant_rnaseq_new_sgRNAs --al coll_trim_cut_contam.fasta --un "$sample_name"_clean.fasta -U coll_trim_cut-noUMI.fasta  > /dev/null

gzip "$sample_name"_clean.fasta

qstat -j $JOB_ID
## end script ##

# submit parallel jobs to wynton
cat samplenames.txt | xargs -n 1 sh -c 'qsub -cwd -pe smp 2 -l mem_free=8G -l h_rt=3:00:00 /wynton/scratch/axu/RP_pre_alignment_withsgRNA_Filter.sh $0'
# took about 1-2 hrs to run on Wynton per sample




# count reads after collapsing (coll_trim_cut.fasta)
###########################
### Count_all_fastas.sh ###
###########################
sample_name=$1

echo "$sample_name"
cd "$sample_name"
grep ">" coll_trim_cut.fasta | wc -l
# end script


cat samplenames.txt | xargs -n 1 sh -c 'bash /wynton/scratch/axu/Count_all_fastas.sh $0'
# sgDDX3_RP_rep1_d0
# 17204584
# sgDDX3_RP_rep1_d14
# 26991185
# sgDDX3_RP_rep1_d3
# 20221545
# sgDDX3_RP_rep2_d0
# 23839591
# sgDDX3_RP_rep2_d14
# 23995459
# sgDDX3_RP_rep2_d3
# 20310045
# sgGal4_RP_rep1_d0
# 31136428
# sgGal4_RP_rep1_d14
# 24455218
# sgGal4_RP_rep1_d3
# 28435513
# sgGal4_RP_rep2_d0
# 34323358
# sgGal4_RP_rep2_d14
# 22554722
# sgGal4_RP_rep2_d3
# 21041862



# count reads after bowtie alignment against repeat RNA (${samplename}_clean.fasta.gz.fasta.gz)
#############################
### Count_all_fasta.gz.sh ###
#############################
sample_name=$1

echo "$sample_name"
cd "$sample_name"
zcat "$sample_name"_clean.fasta.gz | grep ">" | wc -l
# end script


cat samplenames.txt | xargs -n 1 sh -c 'bash Count_all_fasta.gz.sh $0'
sgDDX3_RP_rep1_d0
10318127
sgDDX3_RP_rep1_d14
16909544
sgDDX3_RP_rep1_d3
11724038
sgDDX3_RP_rep2_d0
15304518
sgDDX3_RP_rep2_d14
14860073
sgDDX3_RP_rep2_d3
12543459
sgGal4_RP_rep1_d0
20448780
sgGal4_RP_rep1_d14
15541169
sgGal4_RP_rep1_d3
17827444
sgGal4_RP_rep2_d0
22707270
sgGal4_RP_rep2_d14
14879425
sgGal4_RP_rep2_d3
12676490

# testing L -15
sgDDX3_RP_rep1_d0
10250710
sgDDX3_RP_rep1_d14
16791252
sgDDX3_RP_rep1_d3
11556530
sgDDX3_RP_rep2_d0
15204856
sgDDX3_RP_rep2_d14
14761581
sgDDX3_RP_rep2_d3
12429505
sgGal4_RP_rep1_d0
20307745
sgGal4_RP_rep1_d14
15434173
sgGal4_RP_rep1_d3
17684592
sgGal4_RP_rep2_d0
22557542
sgGal4_RP_rep2_d14
14789761
sgGal4_RP_rep2_d3
12577019

# testing L-10
sgDDX3_RP_rep1_d0
10249079
sgDDX3_RP_rep1_d14
16783265
sgDDX3_RP_rep1_d3
11551040
sgDDX3_RP_rep2_d0
15194156
sgDDX3_RP_rep2_d14
14759222
sgDDX3_RP_rep2_d3
12413328
sgGal4_RP_rep1_d0
20293715
sgGal4_RP_rep1_d14
15433929
sgGal4_RP_rep1_d3
17680174
sgGal4_RP_rep2_d0
22548643
sgGal4_RP_rep2_d14
14785611
sgGal4_RP_rep2_d3
12574159



###################################
### generate_STAR_indices_29.sh ###
###################################
#!/bin/env bash
module load CBI star/2.7.5a

set -x
STAR --runThreadN $NSLOTS \
--runMode genomeGenerate \
--genomeDir /wynton/home/floorlab/axu/references/STAR_indices_29 \
--genomeFastaFiles /wynton/home/floorlab/axu/references/GRCh38.primary_assembly_noscaffold.genome.fa \
--sjdbGTFfile /wynton/home/floorlab/axu/references/gencode.v25.annotation.gtf \
--sjdbOverhang 29 # set as read length - 1

qstat -j $JOB_ID
# end

# 3 hrs 6 core 8G ram



#######################
### RP_alignment.sh ###
#######################

# Uses STAR to make alignments
#!/bin/env bash
set -x

module load CBI star/2.7.5a
module load CBI samtools/1.10
sample_name=$1
cd $sample_name

# lorenzo's paramters
STAR --genomeDir /wynton/home/floorlab/axu/references/STAR_indices_29/ \
--readFilesCommand zcat \
--readFilesIn "$sample_name"_clean.fasta.gz \
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
--alignEndsType Local

samtools index -@ $NSLOTS "$sample_name"_Aligned.sortedByCoord.out.bam

qstat -j $JOB_ID
# end <1hr

# submit parallel jobs to wynton
cat samplenames.txt | xargs -n 1 sh -c 'qsub -cwd -pe smp 2 -l mem_free=4G -l h_rt=3:00:00 /wynton/scratch/axu/RP_alignment.sh $0'


#misc organization

mkdir STAR_logs
find . -type f -name "*.final.out" -exec cp {} STAR_logs \;
scp -r axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu/AXSNF07_concat/STAR_logs/ .

mkdir all_bams
find . -type f -name "*.bam*" -exec mv {} all_bams \;
scp -r axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu/AXSNF07_concat/all_bams/ .

# can run multiqc this on STAR logs for summarizing / generating html reports
# run in folder containint star logs
pip3 install multiqc
/Users/albert/Library/Python/3.8/bin/multiqc .


# cleanup qsub job submissions
find . -type f -name "*.sh.*" -exec mv {} qsub_logs \;


###########################
### Count_all_fastqs.sh ###
###########################
#!/bin/bash
# NOTE : Quote it else use array to avoid problems #
FILES="/bfd/axu/raw_sequencing/AXSNF08_concat/*"
for f in $FILES
do
  echo "Processing $f file..."
  echo "$f"
  # take action on each file. $f store current file name
  zcat "$f" | echo $((`wc -l`/4))
done

###########################
### Count_all_fastqs.sh ###
###########################
#!/bin/bash
# NOTE : Quote it else use array to avoid problems #
FILES="/wynton/scratch/axu/AXSNF07_concat/*"
for f in $FILES
do
  echo "Processing $f file..."
  echo "$f"
  # take action on each file. $f store current file name
  zcat "$f" | echo $((`wc -l`/4))
done

# Processing /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep1_d0.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep1_d0.fastq.gz
# 33283577
# Processing /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep1_d14.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep1_d14.fastq.gz
# 54774530
# Processing /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep1_d3.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep1_d3.fastq.gz
# 45128245
# Processing /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep2_d0.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep2_d0.fastq.gz
# 47468056
# Processing /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep2_d14.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep2_d14.fastq.gz
# 52071546
# Processing /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep2_d3.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgDDX3_RP_rep2_d3.fastq.gz
# 40649002
# Processing /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep1_d0.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep1_d0.fastq.gz
# 67354533
# Processing /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep1_d14.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep1_d14.fastq.gz
# 50508774
# Processing /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep1_d3.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep1_d3.fastq.gz
# 65138611
# Processing /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep2_d0.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep2_d0.fastq.gz
# 77945595
# Processing /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep2_d14.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep2_d14.fastq.gz
# 46201276
# Processing /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep2_d3.fastq.gz file...
# /wynton/scratch/axu/AXSNF07_concat/sgGal4_RP_rep2_d3.fastq.gz
# 44405299







########################
### RNA_alignment.sh ###
########################

# Uses STAR to make alignments
#!/bin/env bash
set -x

module load CBI star/2.7.5a
module load CBI samtools/1.10
sample_name=$1
cd "$sample_name".fastq

# lorenzo's paramters
STAR --genomeDir /wynton/home/floorlab/axu/references/STAR_indices_49/ \
--readFilesCommand zcat \
--readFilesIn "$sample_name".fastq.gz \
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
--alignEndsType Local

samtools index -@ $NSLOTS "$sample_name"_Aligned.sortedByCoord.out.bam

qstat -j $JOB_ID
# end <1hr

# submit parallel jobs to wynton
cat samplenames.txt | xargs -n 1 sh -c 'qsub -cwd -pe smp 4 -l mem_free=4G -l h_rt=24:00:00 RNA_alignment.sh $0'

mkdir STAR_logs
find . -type f -name "*.final.out" -exec cp {} STAR_logs \;

mkdir all_bams_RNA
find . -type f -name "*.bam*" -exec cp {} all_bams_RNA \;

#download files
scp -r axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu/AXSNF06/STAR_logs/ .
scp -r axu@dt2.wynton.ucsf.edu:/wynton/scratch/axu/AXSNF06/all_bams_RNA/ .
