#!/usr/bin/env python
# coding: utf-8

# ## Write final GS MAP / CIs from top and bottom to csv
# Author: Robert Ietswaart  
# Date: 20240118 
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4
# 
# Source: `Timescale_Bayes_20240118_i_prerun.py`  
# For RNA flow Bayesian rate estimation: generate a genes_w_nucdeg_reruns.tsv file that need reruns


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
    k_bound_lo = 1e-4 #1e-4 unit: min^-1: 1 per 7 days
    k_bound_hi = 1e4 #unit: min^-1: 1 per 6 ms, if too restrictive: increase to 1e6
    k_default = [k_bound_lo, k_bound_lo, k_bound_lo, k_bound_hi]
    
    parser = argparse.ArgumentParser(
        description='Write to file the list of genes with GS output for rates.')
    parser.add_argument('--organism',type=str, default='h',
                        help='row index of starting gene for batch.')

    args = parser.parse_args()

    o = args.organism

    inpath = dict()
    inpath['h'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240119_K562')
    inpath['m'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240119_3T3')
    
    outpath = dict()
    outpath['h'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240120_K562')
    outpath['m'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240120_3T3')

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('Bayes_i_prerun')
    log_file = os.path.join(outpath[o],'LogErr', 'Bayes20240120_iii_prerun.log')
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

    organisms = ['m', 'h']
    org_map = {'m': 'mouse', 'h': 'human'}
    fracs = ['chr', 'nuc', 'cyto', 'poly', 'tot']
    fracs_model = {'chr': ['chr_fit'],  
                   'nuc': ['nuc_fit', 
                           'nuc_fit_from_chr_nucres', 
                           'nuc_fit_from_chr_nucdeg'],
                   'cyto': ['cyto_fit_from_nuc', 
                            'cyto_pred_from_chr_nucres', 
                            'cyto_pred_from_chr_nucdeg'],
                   'poly': ['poly_fit_from_nuc'],
                   'tot': ['tot_fit', 'tot_pred_from_nuc', 
                           'tot_pred_from_chr_nucres', 
                           'tot_pred_from_chr_nucdeg']}
    reps = ['G', 'H', 'R', 'S', 'T', 'U']
    red_reps = ['G_R', 'H_S', 'T', 'U']
    org_reps = {'m': ['G','H','R','S'], 'h': ['T', 'U']}
    org_red_reps = {'m': ['G_R','H_S'], 'h': ['T', 'U']}
    red_r = {'G': 'G_R', 'H': 'H_S', 'R': 'G_R', 'S': 'H_S', 'T': 'T', 'U': 'U'}
    red2reps = {'G_R': ['R','G'], 'H_S': ['S','H'], 'T':['T'], 'U':['U']}
    
    OUT_TYPES = ['.MAP', '.0.025.quantile', '.0.975.quantile'] #'.Mean', 
    rt = 'k_'
    POST_PARA_TYPES = ['var_scale', 'underflow']
    
    Timescales = ['T_chr',
                  'T_nuc',
                  'T_cyto',
                  'T_whole_cell',
                  'T_nucdeg',
                  'T_nucexp_from_nucdeg',              
                  'T_chr_release_from_nucdeg']
    
    filename = 'Bayes_Rates_20240119_' + org_map[o] + '.tsv'
    B = pd.read_csv(os.path.join(inpath[o], filename), sep='\t') 
    logger.info(org_map[o])
    logger.info('Load rates that needs inspection of nucdeg rates file \n %s' % filename)
    
    out = []
    for idx in B.index:
        ensid = B['Gene'][idx]
        RERUN= False
        for rr in org_red_reps[o]:
            if not RERUN:
                for ts in Timescales:
                    if ts in {'T_chr', 'T_nuc', 'T_cyto', 'T_whole_cell'}:
                        rate_value = B[rr+'.'+ts.replace('T_',rt)+OUT_TYPES[0]][idx]
                        if np.isnan(rate_value):
                            break
                    if ts in {'T_nucdeg', 'T_nucexp_from_nucdeg', 'T_chr_release_from_nucdeg'}:
                        for i, suf in enumerate(OUT_TYPES):
                            rate_value = B[rr+'.'+ts.replace('T_',rt)+suf][idx]
                            if np.isnan(rate_value):
                                RERUN = True
                                logger.info('%s %s %s RERUN: %s rate is NA' % (ensid, ts, suf, rr))
                                break
                            else:
                                if rate_value != k_default[i]:
                                    break
                                if suf == OUT_TYPES[-1]:#all the other OUT_TYPES are default as well
                                    logger.info('%s %s RERUN: %s rates = default domain, so no information' % (ensid, ts, rr))
                                    RERUN = True
                        if RERUN:
                            break
        if RERUN:
            out.append(ensid)

    out = pd.DataFrame(out)
    filename = 'genes_w_nucdeg_reruns.tsv'
    out.to_csv(os.path.join(outpath[o], filename), sep='\t', index=False, header=False)

    logger.info('end')

    
if __name__ == '__main__':
    main()
