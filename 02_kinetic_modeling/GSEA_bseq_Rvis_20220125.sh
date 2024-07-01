#!/bin/bash
#Robert Ietswaart 
#batch job to visualize GSEA on compartmental Timescales

date=20220127 #first run: 20220125, rerun1:20220126 rerun2:20220127  

program=Rscript
baseDir=/Users/ri23/Documents/HMS/Brendan/GSEA
scDir=/Users/ri23/Documents/HMS/Code
Treatment="chr nuc cyto poly wc"
GENESETS="h.all c5.all c3.all c2.all" #https://www.gsea-msigdb.org/gsea/downloads.jsp
SUFFIX="all filter"

cd $baseDir
mkdir -p ${baseDir}/LogErr

for suf in $SUFFIX
do
    for tr in $Treatment
    do
        for gs in $GENESETS
        do
            project_name=Timescales_GSEA_h_${suf}_${tr}_z_${gs} 
            #first run: Timescales_GSEA_h_${suf}_${tr}_${gs}
            #rerun1: Timescales_GSEA_h_${suf}_${tr}_log_${gs}
            
            echo $project_name
          
            cd ${baseDir}       
          
            ln -s ./${project_name}.GseaPreranked.*/gsea_report_for_na_pos_*.tsv ${project_name}_pos_results.tsv
            ln -s ./${project_name}.GseaPreranked.*/gsea_report_for_na_neg_*.tsv ${project_name}_neg_results.tsv

            Rscript ${scDir}/GSEA_draw_enrichment_figures_bseq20220125.R ${project_name}_pos_results.tsv GSEA 100 80 10 \
            >> ${baseDir}/LogErr/${project_name}_${date}.log    
            #40 30 10: max sig num from timescales in plot (40), max.term.length (30), and point size of text (10)
            #if fewer are significant: then fewer shown
        done 
    done
done
  

