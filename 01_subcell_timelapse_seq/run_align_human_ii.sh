#!/bin/bash

#SBATCH -c 1
#SBATCH -t 60
#SBATCH -p priority
#SBATCH -o %j.out
#SBATCH -e %j.err

# load required modules
module load gcc/6.2.0
module load star/2.7.0a
module load samtools/1.9
module load java/jdk-1.8u112
module load bedtools/2.27.1

# current directory
dir=`pwd`

# set cutadapt dir
cutadapt_dir='/n/groups/churchman/bms36/2021-05-16_T_U/cutadapt'

# list lib names
#files='G1_poly G2_poly G3_poly G4_poly G5_poly H1_poly H2_poly H3_poly H4_poly H5_poly G1_tot G2_tot G3_tot G4_tot G5_tot H1_tot H2_tot H3_tot H4_tot H5_tot'
#files='G1_nuc G2_nuc G3_nuc G4_nuc G5_nuc H1_nuc H2_nuc H3_nuc H4_nuc H5_nuc G1_cyto G2_cyto G3_cyto G4_cyto G5_cyto H1_cyto H2_cyto H3_cyto H4_cyto H5_cyto'
#files='T1_nuc T2_nuc T3_nuc T4_nuc T5_nuc U1_nuc U2_nuc U3_nuc U4_nuc U5_nuc T1_tot T2_tot T3_tot T4_tot T5_tot U1_tot U2_tot U3_tot U4_tot U5_tot'
files='T1_chr T2_chr T3_chr T4_chr T5_chr U1_chr U2_chr U3_chr U4_chr U5_chr T1_cyto T2_cyto T3_cyto T4_cyto T5_cyto U1_cyto U2_cyto U3_cyto U4_cyto U5_cyto T1_poly T2_poly T3_poly T4_poly T5_poly U1_poly U2_poly U3_poly U4_poly U5_poly'
#files='DDX3X_0_tot_rep1 DDX3X_0_tot_rep2 DDX3X_60_chr_rep1 DDX3X_60_chr_rep2 DDX3X_60_cyto_rep1 DDX3X_60_cyto_rep2 DDX3X_60_nuc_rep1 DDX3X_60_nuc_rep2 DDX3X_60_tot_rep1 DDX3X_60_tot_rep2 PABPC4_0_tot_rep1 PABPC4_0_tot_rep2 PABPC4_60_chr_rep1 PABPC4_60_chr_rep2 PABPC4_60_cyto_rep1 PABPC4_60_cyto_rep2 PABPC4_60_nuc_rep1 PABPC4_60_nuc_rep2 PABPC4_60_tot_rep1 PABPC4_60_tot_rep2 scramble_0_tot_rep1 scramble_0_tot_rep2 scramble_60_chr_rep1 scramble_60_chr_rep2 scramble_60_cyto_rep1 scramble_60_cyto_rep2 scramble_60_nuc_rep1 scramble_60_nuc_rep2 scramble_60_tot_rep1 scramble_60_tot_rep2'
#files='DIS3_0_tot_rep1 DIS3_0_tot_rep2 DIS3_60_nuc_rep1 DIS3_60_nuc_rep2 DIS3_60_tot_rep1 DIS3_60_tot_rep2 EXOSC10_0_tot_rep1 EXOSC10_0_tot_rep2 EXOSC10_60_nuc_rep1 EXOSC10_60_nuc_rep2 EXOSC10_60_tot_rep1 EXOSC10_60_tot_rep2 PABPN1_0_tot_rep1 PABPN1_0_tot_rep2 PABPN1_60_nuc_rep1 PABPN1_60_nuc_rep2 PABPN1_60_tot_rep1 PABPN1_60_tot_rep2 ZFC3H1_0_tot_rep1 ZFC3H1_0_tot_rep2 ZFC3H1_60_nuc_rep1 ZFC3H1_60_nuc_rep2 ZFC3H1_60_tot_rep1 ZFC3H1_60_tot_rep2 scramble_0_tot_rep1 scramble_0_tot_rep2 scramble_60_nuc_rep1 scramble_60_nuc_rep2 scramble_60_tot_rep1 scramble_60_tot_rep2'
#files='DIS3_0_tot_rep1'

for x in $files

do

cd /n/groups/churchman/bms36/2021-05-16_T_U/STAR_2023/${x}
	
		OUT=align_${x}.sh 		#THIS is the name of the script to be written. Include var x.
		echo "#!/bin/bash" >> $OUT
		echo "#SBATCH -c 4" >> $OUT
		echo "#SBATCH -t 700" >> $OUT
		echo "#SBATCH -p short" >> $OUT 
		echo "#SBATCH --mem=50G" >> $OUT
		echo "#SBATCH -o %j.out" >> $OUT
		echo "#SBATCH -e %j.err" >> $OUT                      # Partition to run in" >> $OUT
		#echo "#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL" >> $OUT
		#echo "#SBATCH --mail-user=x@g.harvard.edu   # Email to which notifications will be sent" >> $OUT
		echo "" >> $OUT				
		echo "module load gcc/6.2.0" >> $OUT
		echo "module load star/2.7.0a" >> $OUT
		echo "module load samtools/1.9" >> $OUT
		echo "module load java/jdk-1.8u112" >> $OUT
		echo "module load bedtools/2.27.1" >> $OUT
		echo "" >> $OUT
		echo "../align_human.sh $x $cutadapt_dir" >> $OUT
		#echo "Rscript --vanilla ~/scripts/o2_grnaSearch2.sh ../../../fastaToSearch/seqtosearch_${i}_split${n}.fa" >> $OUT
		#echo "date '+%A %W %Y %X'" >> $OUT
		echo "" >> $OUT
		chmod u+x $OUT						# MUST change permission, otherwise will not submit
		sbatch $OUT							# Submit the script 
	
done
