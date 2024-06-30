#!/bin/bash
#Robert Ietswaart 
#sbatch jobs for adapter trimming with cutadapt on Nuc Wash subcell RNAseq fastq files
#For subcell timelapse-seq see cutadapt.sh instead
#Run after bcl2fastq_20230912.sh

date=20231002
scDir=/n/groups/churchman/ri23/code/bseq/
cd ${scDir}
suffix="" #default: "" #alternative (not used): _barcode-mismatches0
parameterFiles=NucWash_CG20230908${suffix}_parameters.in
echo $parameterFiles 

baseDir=`grep "Base Directory"            $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
samples=`grep "Sample Names"           $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
indexDir=`grep "Index Directory"            $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
initDir=`grep "Initial Files Directory"     $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
notifEm=`grep "Notification Email"       $parameterFiles | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`

program1=cutadapt

# load modules
module load gcc/6.2.0
module load python/3.7.4
module load ${program1}/2.5

adapter_R1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
adapter_R2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

nthread=4

FILES=`echo -e $samples`
for f in $FILES
do
    for lane in $(seq 1 4) 
    do 
 
        # cutadapt to trim adapters, 
        # -e .2 : allow 20% error rate since read quality is low at 3' end, r
        # -u 3 : trim 3 nt from read1 5' as they arise from template-switching oligo (see SMARTer Stranded Total RNA Sample Prep Kit - HI Mammalian User Manual)
        # -a -A : trim read1 and read2 3' adapter sequences 
        #-m 20 : require minimum read length of 20
        
        # Second cutadapt: Quality trim at both ends, but adding this quality filter to the above command causes quality trimming before adaptor trimming.
        # -q 20,0 : filter for quality > 20 at read1 and > 0 read2 
        #next-seq-trim=20 :   --nextseq-trim 3'CUTOFF NextSeq-specific quality trimming (each read). Trims also dark cycles appearing as high-quality G bases.
#         -q [5'CUTOFF,]3'CUTOFF, --quality-cutoff [5'CUTOFF,]3'CUTOFF
#                         Trim low-quality bases from 5' and/or 3' ends of each
#                         read before adapter removal. Applied to both reads if
#                         data is paired. 5' end is trimmed with the first cutoff, the 3'
#                         end with the second.

        # Remove intermediates        
        sbatch -n ${nthread} --mem=16G -t 0-00:20:00 --job-name=${program1}${suffix}_${f} \
        -o ${baseDir}/LogErr/${program1}${suffix}_${f}${lane}.log -e ${baseDir}/LogErr/${program1}${suffix}_${f}${lane}.err \
        -p short --mail-user=${notifEm} --mail-type=ALL \
        --wrap="${program1} -a ${adapter_R1} -A ${adapter_R2} -e .2 -m 20 -u 3 -j ${nthread} \
        -o ${initDir}${f}_L00${lane}_R1_001_temp.fastq.gz -p ${initDir}${f}_L00${lane}_R2_001_temp.fastq.gz ${initDir}${f}_L00${lane}_R1_001.fastq.gz ${initDir}${f}_L00${lane}_R2_001.fastq.gz; \
        ${program1} -a ${adapter_R1} -A ${adapter_R2} -q 20,0 --nextseq-trim=20 -m 20 -j ${nthread} \
        -o ${initDir}${f}_L00${lane}_R1_001_cutadapt.fastq.gz -p ${initDir}${f}_L00${lane}_R2_001_cutadapt.fastq.gz ${initDir}${f}_L00${lane}_R1_001_temp.fastq.gz ${initDir}${f}_L00${lane}_R2_001_temp.fastq.gz ; \
        rm ${initDir}${f}_L00${lane}_R1_001_temp.fastq.gz; rm ${initDir}${f}_L00${lane}_R2_001_temp.fastq.gz 
        "  
        
    done
done

echo 'end of cutadapt script'

