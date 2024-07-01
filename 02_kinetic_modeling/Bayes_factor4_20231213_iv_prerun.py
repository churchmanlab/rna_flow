#!/usr/bin/env python
# coding: utf-8

# ## Write final GS MAP / CIs from top and bottom to csv
# Author: Robert Ietswaart  
# Date: 20240106 
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4
# 
# Source: `Timescale_Bayes_20231201_iv_prerun.py`  
# For RNA flow Bayes Factor 4 compartment method: make log_iii file with specs for a run iii for genes that have TIMEOUT in run ii 


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

    outpath = dict()
    outpath['h'] = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor4_20231213_K562')
    outpath['m'] = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor4_20231213_3T3')

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('BayesFactor4_iv_prerun')
    log_file = os.path.join(outpath[o],'LogErr', 'BayesFactor4_20231213_iv_prerun.log')
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

    org_map = {'m': 'mouse', 'h': 'human'}
    
    filename = 'BayesFactor20231213_' + org_map[o] + '_iii.log'
    log_iii = pd.read_csv(os.path.join(outpath[o], 'LogErr', filename), sep='\t') 
    logger.info(org_map[o])
    logger.info('\n %s' % log_iii)
    
    log_iv = dict()
    log_iv['suf_iii'] = []
    log_iv['n_done'] = []
    log_iv['start_id'] = []
    log_iv['n_todo'] = []
    log_iv['suf_iv'] = []
    suffix_iv = max(log_iii['suf_iii']) + 1
    for idx in log_iii.index:
        i = log_iii['suf_iii'][idx]
        n_todo_rerun = log_iii['n_todo'][idx]
        start_id_rerun = log_iii['start_id'][idx]
         
        filename = 'Bayes_factor_20231213_'+org_map[o]+'_batch'+str(i)+'.tsv'
        if os.path.exists(os.path.join(outpath[o], filename)):
            df = pd.read_csv(os.path.join(outpath[o], filename), sep='\t')
            _len = len(df)  
        else:
            _len = 0
    
        if _len < n_todo_rerun:
            if _len == 0:#prep 1 for medium queue
                log_iv['suf_iii'].append(i)
                log_iv['suf_iv'].append(suffix_iv)
                log_iv['n_done'].append(_len)
                log_iv['start_id'].append(start_id_rerun + _len)
                log_iv['n_todo'].append(1)
                logger.info('%s %d %d %d' % (o, i, n_todo_rerun, _len))
                suffix_iv += 1
                _len = 1
                
            if _len < n_todo_rerun:#need for rerun2                
                #prep rest for short
                log_iv['suf_iii'].append(i)
                log_iv['suf_iv'].append(suffix_iv)
                log_iv['n_done'].append(_len)
                log_iv['start_id'].append(start_id_rerun + _len)
                log_iv['n_todo'].append(n_todo_rerun - _len)
                suffix_iv += 1
                logger.info('%s %d %d %d' % (o, i, n_todo_rerun, _len))
            
        else:#run iii finished completely. no need for run iv of batch
            logger.info('%s rerun finished: %d %d %d' % (o, i, n_todo_rerun, _len))

    log_iv = pd.DataFrame(log_iv)
    filename = 'BayesFactor20231213_' + org_map[o] + '_iv.log'
    log_iv.to_csv(os.path.join(outpath[o], 'LogErr', filename), sep='\t', index=False)

    logger.info('end')

    
if __name__ == '__main__':
    main()
