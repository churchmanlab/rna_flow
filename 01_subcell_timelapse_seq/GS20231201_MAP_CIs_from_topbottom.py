#!/usr/bin/env python
# coding: utf-8

# ## Write final GS MAP / CIs from top and bottom to csv
# Author: Robert Ietswaart  
# Date: 20231201  
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4
# 
# Source: `GS_20220118_MAP_CIs_from_topbottom.ipynb`  
# For RNA flow Subcellular Timelapse seq. 


import os
import re
import copy
import numpy as np
import pandas as pd
import logging
import argparse
import fit

from __init__ import default_logger_format, default_date_format


def main():
    np.random.seed(12345)

    parser = argparse.ArgumentParser(
        description='Write final GS MAP / CIs from top and bottom to csv.')
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
    outpath['h'] = os.path.join('/n','groups','churchman','ri23','bseq','GS20231201_K562')
    outpath['m'] = os.path.join('/n','groups','churchman','ri23','bseq','GS20231201_3T3')

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('GS_final')
    log_file = os.path.join(outpath[o],'LogErr', 'GS20231201_MAP_CIs.log')
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
    TB = dict()         #TopBottom genes
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

                    for suffix in SUFFICES[tc]:
                        path_tb = os.path.join(path[o][fr],'STAR_2023',r+'1_'+fr)
                        filename_tb = r + '1_' + fr + '_' + tc + suffix + '.MAPs.txt'
                        if tc == 'top1000':
                            tb_key = red_r[r] + fr + tc
                        if tc == 'bottom500':
                            tb_key = red_r[r] + fr + tc + suffix
                        TB[tb_key]= pd.read_csv(os.path.join(path_tb, filename_tb), 
                                                sep='\t', header=None, 
                                                names=['ENS_ID', 'Symbol', 'MAP'])

    ## Process GS top and bottom to get final MAP / CIs
    for fr in fracs: 
        for r in org_reps[o]:
            rr = red_r[r]
            tc0 = 'top1000' 
            if r+background_id[r]+' MAP' in GS[rr+fr+tc0].columns:
                logger.info('%s, %s, %s' % (o, r, fr))
                N_genes = len(GS[rr+fr+tc0])
                GS[rr+fr+'final'] = dict()
                for col in ['Gene', 'Symbol']:
                    GS[rr+fr+'final'][col] = []
                for t in time_id:
                    for ot in OUT_TYPES:
                        col = r+t+ot
                        GS[rr+fr+'final'][col] = []

                for g_idx in GS[rr+fr+tc0].index:
                    if (g_idx % 100) == 0:
                        logger.info('%d / %d' % (g_idx, N_genes))
                    gene = GS[rr+fr+tc0]['Gene'][g_idx]
                    GS[rr+fr+'final']['Gene'].append(gene)
                    GS[rr+fr+'final']['Symbol'].append(GS[rr+fr+tc0]['Symbol'][g_idx])

                    #get NTR timeseries from top1000 and bottom500
                    TC_TYPES_gene = fit.get_tc_types_for_gene_jit(gene, rr, fr, GS.keys(), TB)

                    NTR = dict()
                    for ot in OUT_TYPES:
                        gs_times = [r+t+ot for t in time_id]                            
                        for i in TC_TYPES_gene:
                            tc = TC_from_jit[i]                      
                            if not ot in NTR.keys():
                                NTR[ot] = np.asarray(GS[rr+fr+tc].loc[g_idx,gs_times], 
                                                     dtype='float64')
                            else:
                                if ot in [' Mean',' MAP']:
                                    NTR[ot] = (NTR[ot] + np.asarray(GS[rr+fr+tc].loc[g_idx, gs_times], 
                                                                    dtype='float64')) / 2
                                elif ot == ' 0.05 quantile':
                                    NTR[ot] = np.minimum(NTR[ot],
                                                         np.asarray(GS[rr+fr+tc].loc[g_idx, gs_times], 
                                                                    dtype='float64'))
                                elif ot == ' 0.95 quantile':
                                    NTR[ot] = np.maximum(NTR[ot],
                                                         np.asarray(GS[rr+fr+tc].loc[g_idx, gs_times], 
                                                                    dtype='float64'))                                    
                        for t in time_id:
                            col = r+t+ot
                            GS[rr+fr+'final'][col].append(NTR[ot][(int(t)-1)])

                GS_final_df = pd.DataFrame(GS[rr+fr+'final'])
                filename = 'GS20231201_' + r + '_' + fr + '.tsv'
                GS_final_df.to_csv(os.path.join(outpath[o], filename), sep='\t',index=False)

    logger.info('end')

    
if __name__ == '__main__':
    main()
