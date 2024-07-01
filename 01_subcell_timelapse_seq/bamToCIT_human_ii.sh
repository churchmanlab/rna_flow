#!/bin/bash

#SBATCH -c 1
#SBATCH -t 60
#SBATCH -p priority
#SBATCH -o %j.out
#SBATCH -e %j.err

# load required modules
module load gcc/6.2.0
module load java/jdk-1.8u112

#cd
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/T-chr

# submit job
OUT=run_T-chr_bamToCIT.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 1" >> $OUT
echo "#SBATCH -t 1440" >> $OUT
echo "#SBATCH -p medium" >> $OUT 
echo "#SBATCH --mem=50G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Bam2CIT -p T-chr_noMT.cit /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T1_chr/T1_chr_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T2_chr/T2_chr_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T3_chr/T3_chr_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T4_chr/T4_chr_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T5_chr/T5_chr_sort.human_noMT.bam" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT	

#cd
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/U-chr

# submit job
OUT=run_U-chr_bamToCIT.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 1" >> $OUT
echo "#SBATCH -t 1440" >> $OUT
echo "#SBATCH -p medium" >> $OUT 
echo "#SBATCH --mem=50G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Bam2CIT -p U-chr_noMT.cit /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U1_chr/U1_chr_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U2_chr/U2_chr_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U3_chr/U3_chr_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U4_chr/U4_chr_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U5_chr/U5_chr_sort.human_noMT.bam" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT	

#cd
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/T-cyto

# submit job
OUT=run_T-cyto_bamToCIT.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 1" >> $OUT
echo "#SBATCH -t 1440" >> $OUT
echo "#SBATCH -p medium" >> $OUT 
echo "#SBATCH --mem=50G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Bam2CIT -p T-cyto_noMT.cit /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T1_cyto/T1_cyto_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T2_cyto/T2_cyto_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T3_cyto/T3_cyto_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T4_cyto/T4_cyto_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T5_cyto/T5_cyto_sort.human_noMT.bam" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT	

#cd
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/U-cyto

# submit job
OUT=run_U-cyto_bamToCIT.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 1" >> $OUT
echo "#SBATCH -t 1440" >> $OUT
echo "#SBATCH -p medium" >> $OUT 
echo "#SBATCH --mem=50G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Bam2CIT -p U-cyto_noMT.cit /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U1_cyto/U1_cyto_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U2_cyto/U2_cyto_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U3_cyto/U3_cyto_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U4_cyto/U4_cyto_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U5_cyto/U5_cyto_sort.human_noMT.bam" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT	

#cd
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/T-poly

# submit job
OUT=run_T-poly_bamToCIT.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 1" >> $OUT
echo "#SBATCH -t 1440" >> $OUT
echo "#SBATCH -p medium" >> $OUT 
echo "#SBATCH --mem=50G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Bam2CIT -p T-poly_noMT.cit /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T1_poly/T1_poly_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T2_poly/T2_poly_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T3_poly/T3_poly_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T4_poly/T4_poly_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/T5_poly/T5_poly_sort.human_noMT.bam" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT	

#cd
cd /n/groups/churchman/bms36/2021-05-16_T_U/GS_2023/U-poly

# submit job
OUT=run_U-poly_bamToCIT.sh 		#THIS is the name of the script to be written. Include var x.
echo "#!/bin/bash" >> $OUT
echo "#SBATCH -c 1" >> $OUT
echo "#SBATCH -t 1440" >> $OUT
echo "#SBATCH -p medium" >> $OUT 
echo "#SBATCH --mem=50G" >> $OUT
echo "#SBATCH -o %j.out" >> $OUT
echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
echo "" >> $OUT				
echo "module load gcc/6.2.0" >> $OUT
echo "module load java/jdk-1.8u112" >> $OUT
echo "" >> $OUT
echo "/n/groups/churchman/bms36/programs/GRAND-SLAM_2.0.5d/gedi -e Bam2CIT -p U-poly_noMT.cit /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U1_poly/U1_poly_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U2_poly/U2_poly_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U3_poly/U3_poly_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U4_poly/U4_poly_sort.human_noMT.bam /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/U5_poly/U5_poly_sort.human_noMT.bam" >> $OUT
echo "" >> $OUT
chmod u+x $OUT						# MUST change permission, otherwise will not submit
sbatch $OUT	
