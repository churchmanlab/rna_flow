#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short


Libs="T1_cyto T2_cyto T3_cyto T4_cyto T5_cyto T1_chr T2_chr T3_chr T4_chr T5_chr T1_poly T2_poly T3_poly T4_poly T5_poly U1_cyto U2_cyto U3_cyto U4_cyto U5_cyto U1_chr U2_chr U3_chr U4_chr U5_chr U1_poly U2_poly U3_poly U4_poly U5_poly"
for lib in $Libs
do
    cp process_n_k_turnover_WT_2023_h.sh $lib
    cd $lib
    sbatch < process_n_k_turnover_WT_2023_h.sh
    cd ..
done
