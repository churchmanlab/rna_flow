#!/bin/bash

#SBATCH -c 1
#SBATCH -t 10
#SBATCH -p priority
#SBATCH -o %j.out
#SBATCH -e %j.err

# load required modules
module load gcc/6.2.0
module load java/jdk-1.8u112
module load R/3.5.1

#	
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/T-chr

# submit job
OUT=GS_T-chr.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 4" >> $OUT
echo "#SBATCH -t 720" >> $OUT
echo "#SBATCH -p short" >> $OUT 
echo "#SBATCH --mem=25G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
#echo "#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL" >> $OUT
#echo "#SBATCH --mail-user=x@g.harvard.edu   # Email to which notifications will be sent" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "module load R/3.5.1" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -no4sUpattern T1 -prefix T_chr -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads T-chr_noMT.cit" >> $OUT
echo "rm T_chr.tsv" >> $OUT
echo "rm T_chr.ext.tsv" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -prefix T_chr -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads T-chr_noMT.cit" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT							# Submit the script 

#	
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/T-cyto

# submit job
OUT=GS_T-cyto.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 4" >> $OUT
echo "#SBATCH -t 720" >> $OUT
echo "#SBATCH -p short" >> $OUT 
echo "#SBATCH --mem=25G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
#echo "#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL" >> $OUT
#echo "#SBATCH --mail-user=x@g.harvard.edu   # Email to which notifications will be sent" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "module load R/3.5.1" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -no4sUpattern T1 -prefix T_cyto -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads T-cyto_noMT.cit" >> $OUT
echo "rm T_cyto.tsv" >> $OUT
echo "rm T_cyto.ext.tsv" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -prefix T_cyto -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads T-cyto_noMT.cit" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT							# Submit the script 

#	
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/T-poly

# submit job
OUT=GS_T-poly.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 4" >> $OUT
echo "#SBATCH -t 720" >> $OUT
echo "#SBATCH -p short" >> $OUT 
echo "#SBATCH --mem=25G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
#echo "#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL" >> $OUT
#echo "#SBATCH --mail-user=x@g.harvard.edu   # Email to which notifications will be sent" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "module load R/3.5.1" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -no4sUpattern T1 -prefix T_poly -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads T-poly_noMT.cit" >> $OUT
echo "rm T_poly.tsv" >> $OUT
echo "rm T_poly.ext.tsv" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -prefix T_poly -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads T-poly_noMT.cit" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT							# Submit the script 

#	
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/U-chr

# submit job
OUT=GS_U-chr.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 4" >> $OUT
echo "#SBATCH -t 720" >> $OUT
echo "#SBATCH -p short" >> $OUT 
echo "#SBATCH --mem=25G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
#echo "#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL" >> $OUT
#echo "#SBATCH --mail-user=x@g.harvard.edu   # Email to which notifications will be sent" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "module load R/3.5.1" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -no4sUpattern U1 -prefix U_chr -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads U-chr_noMT.cit" >> $OUT
echo "rm U_chr.tsv" >> $OUT
echo "rm U_chr.ext.tsv" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -prefix U_chr -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads U-chr_noMT.cit" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT							# Submit the script 

#	
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/U-cyto

# submit job
OUT=GS_U-cyto.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 4" >> $OUT
echo "#SBATCH -t 720" >> $OUT
echo "#SBATCH -p short" >> $OUT 
echo "#SBATCH --mem=25G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
#echo "#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL" >> $OUT
#echo "#SBATCH --mail-user=x@g.harvard.edu   # Email to which notifications will be sent" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "module load R/3.5.1" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -no4sUpattern U1 -prefix U_cyto -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads U-cyto_noMT.cit" >> $OUT
echo "rm U_cyto.tsv" >> $OUT
echo "rm U_cyto.ext.tsv" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -prefix U_cyto -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads U-cyto_noMT.cit" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT							# Submit the script 

#	
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/U-poly

# submit job
OUT=GS_U-poly.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 4" >> $OUT
echo "#SBATCH -t 720" >> $OUT
echo "#SBATCH -p short" >> $OUT 
echo "#SBATCH --mem=25G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
#echo "#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL" >> $OUT
#echo "#SBATCH --mail-user=x@g.harvard.edu   # Email to which notifications will be sent" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "module load R/3.5.1" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -no4sUpattern U1 -prefix U_poly -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads U-poly_noMT.cit" >> $OUT
echo "rm U_poly.tsv" >> $OUT
echo "rm U_poly.ext.tsv" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Slam -allGenes -full -genomic /home/mc348/.gedi/genomic/K562_dm6.oml -prefix U_poly -mode All -snpConv 0.2 -introns -nthreads 4 -strandness Sense -trim5p 10 -trim3p 5 -progress -plot -D -reads U-poly_noMT.cit" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT							# Submit the script 

