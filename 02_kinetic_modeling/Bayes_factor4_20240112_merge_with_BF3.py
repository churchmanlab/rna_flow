#!/usr/bin/env python
# coding: utf-8

# ## merge batches of Bayes Factor 3 compartment method
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
        description='Merge Bayes Factor 4 and 3 compartment method output into single file.')
    
    args = parser.parse_args()

    outpath = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor20240112')
    
    path = dict()
    path['h'] = dict()
    path['m'] = dict()
    
    path['h']['3'] = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor3_20240110_K562')
    path['m']['3'] = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor3_20240110_3T3')
    path['h']['4'] = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor4_20231213_K562')
    path['m']['4'] = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor4_20231213_3T3')
    
    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('BayesFactor43_merge')
    log_file = os.path.join(outpath,'LogErr', 'Bayes_factor4_20240112_merge_with_BF3_py.log')
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)
    
    organisms = ['m', 'h']
    org_map = {'m': 'mouse', 'h': 'human'}
    org_red_reps = {'m': ['G_R','H_S'], 'h': ['T', 'U']}
 
    T_bf = 100

    fit = dict()
    nucdeg = dict()
    BF = dict()
    for o in organisms:
        logger.info('\n %s' % org_map[o])
        fit[o] = dict()
        filename = 'Bayes_factor4_20231213_'+org_map[o]+'.tsv'
        fit[o]['4'] = pd.read_csv(os.path.join(path[o]['4'], filename), sep='\t')
        
        filename = 'Bayes_factor3_20240110_'+org_map[o]+'.tsv'
        fit[o]['3'] = pd.read_csv(os.path.join(path[o]['3'], filename), sep='\t')

        nucdeg[o] = dict()
        filename = 'Nucdeg4_' + org_map[o] + '.tsv'
        nucdeg[o]['4'] = pd.read_csv(os.path.join(path[o]['4'], filename), sep='\t')
                                     
        filename = 'Nucdeg3_' + org_map[o] + '.tsv'
        nucdeg[o]['3'] = pd.read_csv(os.path.join(path[o]['3'], filename), sep='\t')

        logger.info('Comparison between 4 and 3 compartment BF')
        nucdeg3not4 = set(nucdeg[o]['3']['Gene']).difference(set(nucdeg[o]['4']['Gene']))
        logger.info('in 3, not in 4: %d' % len(nucdeg3not4))
        logger.info(nucdeg3not4)
        nucdeg4not3 = set(nucdeg[o]['4']['Gene']).difference(set(nucdeg[o]['3']['Gene']))
        logger.info('in 4, not in 3: %s' % len(nucdeg4not3))
        logger.info(nucdeg4not3)


        logger.info('Generate merge of 4 and 3 compartment BF file')
        BF[o] = fit[o]['4'].merge(fit[o]['3'], on='Gene',how='outer',sort=True,suffixes=['_4','_3'])
        BF[o].insert(loc = 1, column = 'Symbol', 
                     value = BF[o]['Symbol_4'].where(~BF[o]['Symbol_4'].isna(),BF[o]['Symbol_3']))
        BF[o].drop(['Symbol_3', 'Symbol_4'], axis=1, inplace=True)
        for i, rr in enumerate(org_red_reps[o]):
            BF[o].rename({rr+'.var_scale_chr': rr+'.var_scale_chr_4'}, axis=1, inplace=True)
            df = BF[o][rr+'.bayes_factor_4'] > T_bf
            if i == 0:
                pund = df.where(~BF[o][rr+'.bayes_factor_4'].isna(),
                                BF[o][rr+'.bayes_factor_3'] > T_bf)
            else: 
                pund = pund & df.where(~BF[o][rr+'.bayes_factor_4'].isna(), BF[o][rr+'.bayes_factor_3'] > T_bf)
        for rr in org_red_reps[o]:
            pund.loc[BF[o][rr+'.bayes_factor_4'].isna() & BF[o][rr+'.bayes_factor_3'].isna()] = np.nan 
        
        
        BF[o].insert(loc=2, column = 'PUND', value = pund)   
        BF[o]['PUND'] = BF[o]['PUND'].astype('boolean')
        n_genes_w_bf = len(BF[o][~BF[o]['PUND'].isna()])
        n_punds = sum(BF[o][~BF[o]['PUND'].isna()]['PUND'])
        nucdeg_frac = n_punds / n_genes_w_bf
        logger.info('%s %d' % (o,len(BF[o])))
        logger.info('%s %d genes with BFs determined' % (o, n_genes_w_bf))
        logger.info('%s PUNDs %d PUNDs, fraction: %f' % (o, n_punds, nucdeg_frac))

        BF[o].sort_values(by='Gene', inplace=True)
        filename = 'Bayes_factor_20240112_'+org_map[o]+'_final.tsv'
        logger.info('write all genes to file %s' % filename)
        BF[o].to_csv(os.path.join(outpath, filename), sep='\t', index=False)
                                     
        nucdeg[o]['43'] = BF[o][BF[o]['PUND']==True].sort_values(by='Gene')
        filename = 'Nucdeg_20240112_'+org_map[o]+'_final.tsv'
        logger.info('write PUNDs to file %s' % filename)
        nucdeg[o]['43'].to_csv(os.path.join(outpath, filename), sep='\t', index=False)
    
    logger.info('end')
    
if __name__ == '__main__':
    main()
