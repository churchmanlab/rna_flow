#!/bin/bash
#Robert Ietswaart 
#batch job to run GSEA (v4.2.2 using command line) for ranked standardized half lives (z-scores)

date=20240204

program=/Users/ri23/Documents/HMS/GSEA_4.2.2/gsea-cli.sh
baseDir=/Users/ri23/Documents/HMS/Brendan/GSEA${date}
scDir=/Users/ri23/Documents/HMS/Code
Timescales="chr chr_release nuc nucexp cyto poly_entry whole_cell"

GENESETS="h.all c5.all c3.all c2.all" #https://www.gsea-msigdb.org/gsea/downloads.jsp
OUT_TYPES=.Mean

cd $baseDir
mkdir -p ${baseDir}/LogErr


for ot in $OUT_TYPES
do
    for tr in $Timescales
    do
        input_file=GSEA_h_all_half_life_${tr}${ot}_z 

        for gs in $GENESETS
        do
            project_name=${input_file}_${gs}
            echo $project_name
            
            $program GSEAPreranked -gmx ${baseDir}/msigdb_v7.5.1_GMTs/${gs}.v7.5.1.symbols.gmt \
            -scoring_scheme weighted -norm meandiv -nperm 1000 -rnk ${input_file}.rnk \
            -rpt_label ${project_name} -create_svgs false -make_sets true -plot_top_x 100 -rnd_seed timestamp \
            -set_max 500 -set_min 15 -zip_report false -out ./ >> ${baseDir}/LogErr/${project_name}_${date}.log

       done 
    done
done
