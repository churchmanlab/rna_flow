#!/bin/bash
#Robert Ietswaart 

date=20240126
baseDir=/n/groups/churchman/bms36/2022-11-04_RBP_KDs_ii/
resourceDir=/n/groups/churchman/ri23/bseq/GS${date}/
# scDir=/n/groups/churchman/ri23/code/

notifEm=robert_ietswaart@hms.harvard.edu

nthread=4
athread=$(($nthread-1))

# load required modules
module load gcc/6.2.0
module load java/jdk-1.8u112
module load bedtools/2.27.1
module load samtools/1.9

mkdir -p ${baseDir}/LogErr

prefixes='DIS3_60_nuc_rep1 DIS3_60_nuc_rep2 DIS3_60_tot_rep1 DIS3_60_tot_rep2 EXOSC10_60_nuc_rep1 EXOSC10_60_nuc_rep2 EXOSC10_60_tot_rep1 EXOSC10_60_tot_rep2 PABPN1_60_nuc_rep1 PABPN1_60_nuc_rep2 PABPN1_60_tot_rep1 PABPN1_60_tot_rep2 ZFC3H1_60_nuc_rep1 ZFC3H1_60_nuc_rep2 ZFC3H1_60_tot_rep1 ZFC3H1_60_tot_rep2 scramble_60_nuc_rep1 scramble_60_nuc_rep2 scramble_60_tot_rep1 scramble_60_tot_rep2'
# prefixes=DIS3_60_nuc_rep1 #code dev
# prefixes=EXOSC10_60_nuc_rep1 #unfinished rerun

read_types='unspliced_junctions introns not_introns' #'unspliced_junctions introns' #

for x in $prefixes 
do 
    target_folder=${baseDir}STAR_2023/${x}/  
    input_bam=${target_folder}${x}_sort.human_noMT.bam
    
    for rt in $read_types
    do
        jobname=GS${date}_KDii_extract_${x}_${rt}
        echo $jobname
        
        output_bam=${target_folder}${rt}/${x}_sort.human_noMT_${rt}.bam
        
        if [ $rt == 'unspliced_junctions' ]
        then
            # Filter out split reads, i.e. containing a N in their CIGAR string, as they are not unspliced
            # Select reads from input.bam NOT fully overlapping with exons, 
            # Further select reads that NOT fully overlap with introns.
            # Filter out reads that do not have a mate in the same bam file.
            
            sbatch -p short -n ${nthread} --mem=8G -t 0-01:00:00 --job-name=${jobname} \
            -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
            --mail-user=${notifEm} --mail-type=ALL \
            --wrap="samtools view -@ $athread -h $input_bam | awk '\$1 ~ /^@/ || \$6 !~ /N/' | samtools view -@ $athread -h -b | bedtools intersect -abam stdin -b ${resourceDir}exons.gtf -f 1.0 -v -wa | bedtools intersect -abam stdin -b ${resourceDir}introns.gtf -f 1.0 -v -wa | samtools view -@ $athread -F 0x8 -h -b | samtools sort -@ $athread -o ${output_bam}; samtools index -@ $athread ${output_bam} ${output_bam}.bai"

        elif [ $rt == 'introns' ]
        then
            
             sbatch -p short -n ${nthread} --mem=8G -t 0-01:00:00 --job-name=${jobname} \
            -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
            --mail-user=${notifEm} --mail-type=ALL \
            --wrap="bedtools intersect -a ${input_bam} -b ${resourceDir}introns_sort.bed -f 1.0 -wa | samtools view -@ $athread -F 0x8 -h -b | samtools sort -@ $athread -o ${output_bam}; samtools index -@ $athread ${output_bam} ${output_bam}.bai"       
        
        elif [ $rt == 'not_introns' ]
        then
        
            sbatch -p short -n ${nthread} --mem=8G -t 0-01:00:00 --job-name=${jobname} \
            -o ${baseDir}/LogErr/${jobname}.log -e ${baseDir}/LogErr/${jobname}.err \
            --mail-user=${notifEm} --mail-type=ALL \
            --wrap="bedtools intersect -a ${input_bam} -b ${resourceDir}introns_sort.bed -f 1.0 -wa -v | samtools view -@ $athread -F 0x8 -h -b | samtools sort -@ $athread -o ${output_bam}; samtools index -@ $athread ${output_bam} ${output_bam}.bai"
            
        fi
    done
done
#-p short -n ${nthread} --mem=8G -t 0-01:00:00
#-p short -n ${nthread} --mem=8G -t 0-03:00:00 #unfinished reruns


#NOT USED
# module load gcc/9.2.0
# module load java/jdk-1.8u112
# module load bedtools/2.30.0
# module load samtools/1.15.1 
