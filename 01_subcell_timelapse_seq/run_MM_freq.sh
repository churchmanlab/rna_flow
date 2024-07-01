#!/bin/bash

#SBATCH -c 1
#SBATCH -t 120
#SBATCH -p priority
#SBATCH --mem=100G
#SBATCH -o %j.out
#SBATCH -e %j.err


# load required modules
module load gcc/6.2.0
module load R/3.5.1

./MM_U_chr.R
./MM_T_cyto.R
./MM_U_cyto.R
./MM_T_poly.R
./MM_U_poly.R
./MM_U_chr_slow.R
./MM_T_cyto_slow.R
./MM_U_cyto_slow.R
./MM_T_poly_slow.R
./MM_U_poly_slow.R
