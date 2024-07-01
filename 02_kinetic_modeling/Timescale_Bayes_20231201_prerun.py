#!/usr/bin/env python
# coding: utf-8

# ## Write final GS MAP / CIs from top and bottom to csv
# Author: Robert Ietswaart  
# Date: 20231201  
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4
# 
# Source: `GS20231201_MAP_CIs_from_topbottom.py`  
# For RNA flow Bayesian rate estimation: generate a list of gene ID for which data is available. 


import os
import re
import copy
import numpy as np
import pandas as pd
import logging
import argparse

from __init__ import default_logger_format, default_date_format


def main():
    np.random.seed(12345)

    parser = argparse.ArgumentParser(
        description='Write to file the list of genes with GS output for rates.')
    parser.add_argument('--organism',type=str, default='h',
                        help='row index of starting gene for batch.')

    args = parser.parse_args()

    o = args.organism

    path = dict() 
    path['h'] = dict()
    path['m'] = dict()
    path['h']['chr'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2021-05-16_T_U')
    path['h']['nuc'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2021-04-07_T_U')
    path['h']['cyto'] = path['h']['chr']
    path['h']['poly'] = path['h']['chr']
    path['h']['tot'] = path['h']['nuc']
    path['m']['chr'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2021-04-05_R_S')
    path['m']['nuc'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2020-10-06_G_H_fract')
    path['m']['cyto'] = path['m']['nuc']
    path['m']['poly'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2020-07-13_G_H')
    path['m']['tot'] = path['m']['poly']

    outpath = dict()
    outpath['h'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240119_K562')#'Bayes20231201_K562')
    outpath['m'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240119_3T3')#'Bayes20231201_3T3')

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('Bayes_prerun')
    log_file = os.path.join(outpath[o],'LogErr', 'Bayes20231201_prerun.log')
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

    fracs = ['chr', 'nuc', 'cyto', 'poly', 'tot']
    reps = ['G','H','R','S','T','U']
    org_reps = {'m': ['G','H','R','S'], 'h': ['T', 'U']}
    red_r = {'G': 'G_R', 'H': 'H_S', 'R': 'G_R', 'S': 'H_S', 'T': 'T', 'U': 'U'}
    red_reps = ['G_R', 'H_S', 'T', 'U']
    org_red_reps = {'m': ['G_R','H_S'], 'h': ['T', 'U']}
    time_id = [str(i) for i in range(1,6)]
    background_id = {r: '1' for r in reps}
    # time_mins = [0, 15, 30, 60, 120]

    ptc_genes = dict() #List with protein-coding genes ENS IDs
    ptc_genes['h'] = os.path.join('/n', 'groups', 'churchman', 'mc348','TimelapseSeq','SeqFiles',
                                  'K562_ensGRCh38_MTmod_ptc_list.txt')
    ptc_genes['m'] = os.path.join('/n', 'groups', 'churchman', 'mc348','TimelapseSeq','SeqFiles',
                                  'NIH3T3_mm10_MTmod_ptc_list.txt')
    PTC = pd.read_csv(os.path.join(ptc_genes[o]))

    TC_TYPES = ['top1000','bottom500']#BM model turnover method to estimate TC
    TC_from_jit = {1: 'top1000', 0: 'bottom500'}
    OUT_TYPES = [' 0.05 quantile', ' Mean', ' MAP', ' 0.95 quantile'] #GS specific
    SUFFICES = {'top1000': ['genes.turnover'], 'bottom500': ['', '_below']}

    GS = dict()         #GRAND-SLAM
    for r in reps:
        for fr in fracs:
            for tc in TC_TYPES:
                path_gs = os.path.join(path[o][fr],'GS_2023',r+'-'+fr,re.sub(r'\d', '', tc))
                filename_gs = r + '_' + fr + '.tsv'
                if os.path.exists(os.path.join(path_gs, filename_gs)):
                    GS[red_r[r]+fr+tc]= pd.read_csv(os.path.join(path_gs, filename_gs), sep='\t')
                    #filter for protein-coding genes
                    GS[red_r[r]+fr+tc] = GS[red_r[r]+fr+tc][GS[red_r[r]+fr+tc]['Gene'].isin(PTC['Gene'])]
                    GS[red_r[r]+fr+tc].sort_values(by='Gene',inplace=True, ignore_index=True)

    
    #determine genes that have sufficient data to estimate timescales 
    #NB: top1000 and bottom500 GS have identical gene ids
    genes_w_rates = dict()
    for r in org_reps[o]:
        for fr in fracs:
            if red_r[r]+fr+'top1000' in GS.keys():
                if fr == 'chr':
                    genes_w_rates[red_r[r]] = set(GS[red_r[r]+fr+'top1000']['Gene'])
                elif fr in ['cyto','poly','nuc','tot']:
                    genes_w_rates[red_r[r]] = genes_w_rates[red_r[r]].union(GS[red_r[r]+fr+'top1000']['Gene'])
                logger.info('%s %s %d' % (red_r[r], fr, len(genes_w_rates[red_r[r]])))
    genes_w_rates = genes_w_rates[org_red_reps[o][0]].union(genes_w_rates[org_red_reps[o][1]])
    genes_w_rates = sorted(list(genes_w_rates))
    logger.info('total number of genes: union of %s: %d' % (org_red_reps[o],len(genes_w_rates)))

    genes_df = pd.Series(genes_w_rates)

    filename = 'genes_w_rates.csv'
    genes_df.to_csv(os.path.join(outpath[o], filename), sep='\t',index=False, header=False)
    

    logger.info('end')

    
if __name__ == '__main__':
    main()
