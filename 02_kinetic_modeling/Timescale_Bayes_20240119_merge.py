#!/usr/bin/env python
# coding: utf-8

# ## Write final GS MAP / CIs from top and bottom to csv
# Author: Robert Ietswaart  
# Date: 20231226 
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4
# 
# Source: `Timescale_Bayes_20231201_rerun.py`  
# For RNA flow Bayesian rate estimation: merge the batch output files into one final rates file 


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
        description='Merge Bayes rates output from batches into single file.')
    parser.add_argument('--organism',type=str, default='h',
                        help='organism: h (human K562) or m (mouse 3T3)')
    parser.add_argument('--nbatch',type=int, default=0,
                        help='number of batches to merge')
    
    args = parser.parse_args()

    o = args.organism

    outpath = dict()
    outpath['h'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240119_K562')
    outpath['m'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240119_3T3')

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('Bayes_merge')
    log_file = os.path.join(outpath[o],'LogErr', 'Bayes20240119_merge.log')
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

    org_map = {'m': 'mouse', 'h': 'human'}
    org_red_reps = {'m': ['G_R','H_S'], 'h': ['T', 'U']}
    RATE_TYPES = ['half_life_','k_']
    POST_PARA_TYPES = ['var_scale', 'underflow']

    Timescales = ['T_chr',
                  'T_nuc',
                  'T_nucexp_from_nucres',
                  'T_nucexp_from_dist',
                  'T_cyto',
                  'T_poly_entry',
                  'T_whole_cell',
                  'T_nucdeg',
                  'T_nucexp_from_nucdeg',              
                  'T_chr_release_from_nucdeg']
        
    logger.info('\n %s' % org_map[o])
    fits_all = pd.DataFrame()
    logger.info('append batch files')
    for i in range(args.nbatch):
        filename = 'Bayes_Rates_20240119_'+org_map[o]+'_batch'+str(i)+'.tsv'
        if os.path.exists(os.path.join(outpath[o], filename)):
            df = pd.read_csv(os.path.join(outpath[o], filename), sep='\t') 
            fits_all = fits_all.append(df, ignore_index=True)                    
        else:
            logger.warning('%s batch %d missing' % (o,i))

    logger.info('delete repetitive and all NA columns')
    rt = 'half_life_'
    drop_cols = []
    for rr in org_red_reps[o]:
        for ts in Timescales:
            for out in POST_PARA_TYPES:
                drop_cols.append(rr+'.'+ts.replace('T_',rt)+'.'+out)
    for rt in RATE_TYPES:
        for rr in org_red_reps[o]:
            for ts in Timescales[7:]:#nucdeg/nucexp_from_nucdeg/chr_release_from_nucdeg
                out = 'Mean' #could not be calculated, so were set equal to MAP: redundant info so delete
                drop_cols.append(rr+'.'+ts.replace('T_',rt)+'.'+out)
    fits_all.drop(columns=drop_cols, axis=1, inplace=True)
    fits_all.dropna(axis=1, how='all', inplace=True)
    fits_all.reset_index(drop=True, inplace=True)
    fits_all = fits_all[~fits_all['Gene'].duplicated()]
    
    filename = 'Bayes_Rates_20240119_'+org_map[o]+'.tsv'
    logger.info('write to file %s' % filename) 
    fits_all.sort_values(by='Gene', inplace=True)
    fits_all.to_csv(os.path.join(outpath[o], filename), sep='\t', index=False)
    logger.info('%s %d' % (org_map[o], len(fits_all)))
    
    logger.info('end')

    
if __name__ == '__main__':
    main()
