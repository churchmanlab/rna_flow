#!/usr/bin/env python
# coding: utf-8

# ## merge batches of Bayes Factor 4 compartment method
# Author: Robert Ietswaart  
# Date: 20240111 
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4
# 
# Source: `Bayes_factor4_20231213_merge.py`
# For RNA flow Bayesian Factor calculation: merge the batch output files 


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
        description='Merge Bayes Factor 4 compartment method output from batches into single file.')
    parser.add_argument('--organism',type=str, default='h',
                        help='organism: h (human K562) or m (mouse 3T3)')
    parser.add_argument('--nbatch',type=int, default=0,
                        help='number of batches to merge')
    
    args = parser.parse_args()

    o = args.organism

    outpath = dict()
    outpath['h'] = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor4_20231213_K562_oflow1')
    outpath['m'] = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor4_20231213_3T3')
        
    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('BayesFactor4_merge')
    log_file = os.path.join(outpath[o],'LogErr', 'BayesFactor4_20231213_merge.log')
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

    org_map = {'m': 'mouse', 'h': 'human'}
    org_red_reps = {'m': ['G_R','H_S'], 'h': ['T', 'U']}
    
    logger.info('\n %s' % org_map[o])
    fits_all = pd.DataFrame()
    logger.info('append batch files')
    for i in range(args.nbatch):
        filename = 'Bayes_factor_20231213_'+org_map[o]+'_batch'+str(i)+'.tsv'
        if os.path.exists(os.path.join(outpath[o], filename)):
            df = pd.read_csv(os.path.join(outpath[o], filename), sep='\t') 
            fits_all = fits_all.append(df, ignore_index=True)                    
        else:
            logger.warning('%s batch %d missing' % (o,i))

    fits_all.reset_index(drop=True, inplace=True)
    fits_all = fits_all[~fits_all['Gene'].duplicated()]
    fits_all.sort_values(by='Gene', inplace=True)
    
    filename = 'Bayes_factor4_20231213_'+org_map[o]+'_controli.tsv'
    logger.info('write to file %s' % filename) 
    fits_all.to_csv(os.path.join(outpath[o], filename), sep='\t', index=False)
    logger.info('%s %d' % (org_map[o], len(fits_all)))
    
    T_bf = 100
    nucdeg = copy.deepcopy(fits_all)
    genes_w_data = copy.deepcopy(fits_all)
    for r in org_red_reps[o]:
        nucdeg = nucdeg[~nucdeg[r+'.bayes_factor'].isna()]
        nucdeg = nucdeg[nucdeg[r+'.bayes_factor'] > T_bf]
        genes_w_data = genes_w_data[~genes_w_data[r+'.bayes_factor'].isna()]
    nucdeg_frac = len(nucdeg) / len(genes_w_data)
    
    filename = 'Nucdeg4_'+org_map[o]+'_controli.tsv'
    nucdeg.to_csv(os.path.join(outpath[o], filename), sep='\t', index=False)
    logger.info('Write to file %s %s: %d / %d = %f ' % (org_map[o], filename, len(nucdeg), len(genes_w_data), nucdeg_frac))
    
    logger.info('end')
    
if __name__ == '__main__':
    main()
