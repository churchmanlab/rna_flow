#!/usr/bin/env python
# coding: utf-8

# ## Machine learning to explain rates
# Author: Robert Ietswaart  
# Date: 20230329  
# License: BSD2.  
# Python v3.7.4
# 
# Source: `ML_20220618_subcell_feat_select1.py`  
# For RNA flow project.
# See also corresponding slurm batch script which calls this .py script: 
# ML_20230329_subcell_feat_select1.sh

import os
import re
import copy
import numpy as np
import pandas as pd
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
import logging
import argparse

from sklearn.model_selection import train_test_split, KFold
from sklearn.ensemble import RandomForestRegressor
from sklearn import linear_model

from __init__ import default_logger_format, default_date_format


def main():

    parser = argparse.ArgumentParser(
        description='Machine learning feature selection round 1 on Bayes rates.')
    parser.add_argument('--rate',type=str, default='whole_cell',help='provide a rate type. Choose from: '
                        'chr, chr_release, nuc, nucexp, cyto, poly_entry, whole_cell')
    parser.add_argument('--feats',type=int, default=1,
                        help='tested feature file(s) numberical index. Choose from:'
                             '1: gene_structure_v2.txt gene_structure_agarwal.txt '
                             '2: gene_sequence_v2.txt gene_sequence_agarwal.txt '
                             '3: Codon_features.tsv.gz 4: kmer_ORF_features.tsv.gz'
                             '5: kmer_5UTR_features.tsv.gz 6: kmer_3UTR_features.tsv.gz'
                             '7: polyA_tail_length_data.txt 8: RNA_modifications_v2.txt '
                             '9: gene_location_v2.txt 10: histone_modifications_v1.txt '
                             '11: microRNA_targets_v2.txt 12: RBP_targets_filter_v2.txt '
                             '13: TF_scRNAseq_v1.csv 14: TF_MSigDB_targets_v2.tsv')

    num2feat_files = {1: 'gene_structure_v2.txt gene_structure_agarwal.txt',
                      2: 'gene_sequence_v2.txt gene_sequence_agarwal.txt',
                      3: 'Codon_features.tsv.gz', 4: 'kmer_ORF_features.tsv.gz',
                      5: 'kmer_5UTR_features.tsv.gz', 6: 'kmer_3UTR_features.tsv.gz',
                      7: 'polyA_tail_length_data.txt', 8: 'RNA_modifications_v2.txt',
                      9: 'gene_location_v2.txt', 10: 'histone_modifications_v1.txt',
                      11: 'microRNA_targets_v2.txt', 12: 'RBP_targets_filter_v2.txt',
                      13: 'TF_scRNAseq_v1.csv', 14: 'TF_MSigDB_targets_v2.tsv'}

    args = parser.parse_args()

    path = os.path.join('/n','groups','churchman','ri23','bseq','ML20230329')
    feature_path = os.path.join('/n','groups','churchman','ri23','bseq','RF20220426','features')

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('ML')
    log_file = os.path.join(path,'LogErr', 'ML_20230329_subcell_select1_k_%s_feat%s_py.log' % (args.rate, str(args.feats)))
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

    
    logger.info('Load rates: dependent variables in model')
    organisms = ['m','h']
    org_map = {'m': 'mouse', 'h': 'human'}
    org_red_reps = {'m': ['G_R','H_S'], 'h': ['T', 'U']}

    k_bound_lo = 1e-4 #unit: min^-1: 1 per 7 days
    k_bound_hi = 1e4 #unit: min^-1: 1 per 6 ms

    RATE_TYPE = ['half_life_','k_']

    rt = RATE_TYPE[1] #ML model dependent variable: rates

    Timescales = ['chr',
                  'chr_release',
                  'nucdeg',
                  'nucexp',
                  'nuc',
                  'cyto',
                  'poly_entry',
                  'whole_cell']
    Timescales = [rt + ts for ts in Timescales]

    OUT_TYPES = ['.Mean', '.MAP', '.0.975.quantile', '.0.025.quantile']

    B = dict()          #Bayes fits file
    K = dict()          #Bayes Factor

    o = 'h'
    # for o in organisms:    
    path_b = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20230128')
    filename_b = 'Bayes_Rates_20230128_'+ org_map[o] + '.tsv'
    path_k = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor20221206')
    filename_k = 'Bayes_factor_20230317_' + org_map[o] + '_final.tsv'

    B[o] = pd.read_csv(os.path.join(path_b, filename_b), sep='\t')
    K[o] = pd.read_csv(os.path.join(path_k, filename_k), sep='\t')

    logger.info('Preprocess rates')
    # - preprocess nucexp / chr_release: depending on PUND no or yes 
    # - clip to domain bounds
    # - log transform
    # - standardize


    T_bf = 100

    C = copy.deepcopy(B)
    # for o in organisms:
    C[o] = C[o].merge(K[o], on='Gene', how='outer', suffixes=('', '_bf'))
    logger.info('%s all genes: %d' % (org_map[o], len(C[o])))

    for ot in OUT_TYPES:
        for rr in org_red_reps[o]:
            ts = rt + 'chr_release'
            C[o][rr+'.'+ts+ot] = copy.deepcopy(C[o][rr+'.'+ts+'_from_nucdeg'+ot].where(
                    C[o]['PUND'], C[o][rr+'.'+rt+'chr'+ot]))
            ts = rt + 'nucdeg'#only for nucdeg genes according to Bayes Factor
            C[o][rr+'.'+ts+ot].where(C[o]['PUND'], np.nan, inplace=True)
            ts = rt + 'nucexp'
            C[o][rr+'.'+ts+ot] = copy.deepcopy(C[o][rr+'.'+ts+'_from_nucdeg'+ot].where(
                C[o]['PUND'], C[o][rr+'.'+ts+'_from_nucres'+ot])) 


            for ts in Timescales:
                ###Clip range of values beyond numerical integration domain bounds                  
                C[o][rr+'.'+ts+ot].where(((C[o][rr+'.'+ts+ot] > k_bound_lo) | (C[o][rr+'.'+ts+ot].isna())), 
                                         k_bound_lo, inplace=True) 
                C[o][rr+'.'+ts+ot].where(((C[o][rr+'.'+ts+ot] < k_bound_hi) | (C[o][rr+'.'+ts+ot].isna())), 
                                         k_bound_hi, inplace=True)  

                #standardize rates (to Z-score) to enhance learning
                rates = np.log(C[o][rr+'.'+ts+ot])
                C[o][rr+'z'+ts+ot] = (rates - rates.mean()) / rates.std()

    
    rand_seed = 42
    shuffle = True
    test_size = 0.1
    logger.info('Master Train/CV vs %s test set split' % test_size)
    train_idx, test_idx = train_test_split(C[o].index, 
                                           test_size = test_size, 
                                           random_state = rand_seed, 
                                           shuffle = shuffle)
    C[o] = C[o].loc[train_idx,:]


    logger.info('Load features')
    
    F_raw = dict()

    feature_files = num2feat_files[args.feats].split()

    for f in feature_files:
        F_raw[f] = pd.read_csv(os.path.join(feature_path, f), sep='\t')

        logger.info('%s %d' % (f, len(F_raw[f])))

        if not (f in {'Codon_features.tsv.gz',#too many features to list
                      'kmer_ORF_features.tsv.gz',
                      'kmer_5UTR_features.tsv.gz',
                      'kmer_3UTR_features.tsv.gz'}):
            logger.info(F_raw[f].columns)
            logger.info(F_raw[f].head())   


    logger.info('Generate feature matrix')
    
    def get_merge_key(cols,f):
        return 'Gene'

    def get_cols_merge(cols):
        cols = list(cols)
        if 'Gene' in cols:
            if 'Symbol' in cols:
                cols.remove('Symbol')
        return cols

    for i, f in enumerate(feature_files):
        cols = F_raw[f].columns
        sm_key = get_merge_key(cols,f)#start or merge key
        cols_merge = get_cols_merge(cols)
        F_raw[f] = F_raw[f][F_raw[f][sm_key].isin(C[o][sm_key])]

        if f in {'gene_structure_v2.txt', 'gene_sequence_agarwal.txt',
                 'gene_sequence_v2.txt', 'gene_structure_agarwal.txt',
                 'Codon_features.tsv.gz', 'kmer_ORF_features.tsv.gz',
                 'kmer_5UTR_features.tsv.gz', 'kmer_3UTR_features.tsv.gz',
                 'polyA_tail_length_data.txt', 'RNA_modifications_v2.txt'}:
            logger.info('%s merge with standardization' % f) 

            #standardize quantitative features (get Z-score) to enhance learning
            feats = copy.deepcopy(F_raw[f][cols_merge])
            feats[cols_merge[1:]] = (feats[cols_merge[1:]] - feats[cols_merge].mean()) / feats[cols_merge].std()

        elif f in {'gene_location_v2.txt', 'histone_modifications_v1.txt',
                   'microRNA_targets_v2.txt', 'RBP_targets_filter_v2.txt',
                   'TF_scRNAseq_v1.csv', 'TF_MSigDB_targets_v2.tsv'}:
            logger.info('%s pivot, one hot encode %s' % (f, cols_merge))#no standardization needed

            feats = F_raw[f][cols_merge].pivot(index=cols_merge[0], 
                                               columns=cols_merge[1], 
                                               values=cols_merge[1:])

            for i, c in enumerate(cols_merge[1:]):
                if i == 0:
                    feats[c] = feats[c].where(feats[c].isna(), 1)
                feats[c] = feats[c].where(~feats[c].isna(), 0)
            feats.columns = ['_'.join(col) for col in feats.columns.values] #from multi-index to index     
            feats.reset_index(inplace=True) #Gene no longer index
            
            feats = C[o][[sm_key]].merge(feats, on=sm_key, how ='left', suffixes=('', f))#feat row now all rates

        else:
            logger.info('fail: %s not recognized' % f)

        if i == 0:   
            F = copy.deepcopy(feats)
            logger.info('%s initial feature mat: #features: %d ,#genes: %d' % (f, len(F.columns)-1, len(F))) 
        else:
            F = F.merge(feats, on=sm_key, how ='outer', suffixes=('', f))
            logger.info('%s past merge feature mat: #features: %d ,#genes: %d' % (f, len(F.columns)-1, len(F)))


    #Missing data imputation: NaN -> 0
    #Correct for one hot encoded: no target gene, since only target genes are included in data files,
    #so genes not in file, i.e. NaN genes, are not target genes by construction
    #Least biased imputation if standardized feature: 0 = mean of variable across data.
    F = F.where(~F.isna(), 0)

    #Filter for genes present in RNA flow data
    F = F[F['Gene'].isin(C[o]['Gene'])]
    n_features = len(F.columns) - 1
    logger.info('Final feature matrix: #features: %d ,#genes: %d' % (n_features, len(F)))

    #Get a mapping for features to column number
    feat2id = dict()
    id2feat = dict()
    for i, c in enumerate(F.columns[2:]):#skip Gene and Symbol
        feat2id[c] = i
        id2feat[i] = c


   
    logger.info('Cross validation and non zero feature selection, for each replicate separately')

    ts = rt + args.rate
    ot = OUT_TYPES[0] #NB: you tested predictive power is higher when predicting .Mean over .MAP. 
    logger.info('LASSO1 CV on rate %s' % ts+ot)

    LS = dict()
    fitted_genes = dict()

    if n_features < 1000:
        L1_alpha = [0.0001 * 10**(i) for i in range(4)]
    else:
        L1_alpha = [0.001 * 10**(i) for i in range(3)]

    n_cores = 1
    cv_folds = 10
    max_iter_ = 2e4 #default 1e3 (take 3mins), 1e6 can take >1.5h 
    logger.info('n CV folds: %d' % cv_folds)

    for rr in org_red_reps[o]:
        logger.info('rep %s' % (rr))

        genes_w_rates = ~C[o][rr+'z'+ts+ot].isna()

        rates = copy.deepcopy(C[o][genes_w_rates][['Gene', rr+'z'+ts+ot]])
        rates = rates[rates['Gene'].isin(F['Gene'])]
        rates.sort_values(by=['Gene'], inplace=True)
        rates.reset_index(inplace=True, drop=True)

        features = copy.deepcopy(F[F['Gene'].isin(C[o][genes_w_rates]['Gene'])])
        features = features[features['Gene'].isin(F['Gene'])]
        features.sort_values(by=['Gene'], inplace=True)
        features.reset_index(inplace=True, drop=True)

        check1 = sum(rates['Gene'] == features['Gene'])
        logger.info('Confirm that rates and features are same order: '
                    '%d = %d' % (len(rates), check1))
        fitted_genes[rr] = rates['Gene']
        rates.drop(['Gene'], axis=1, inplace=True)
        features.drop(['Gene'], axis=1, inplace=True) 

        #Lasso for feature selection
        LS['r2_train'] = dict()
        LS['r2_cv'] = dict()
        LS['f'] = dict()

        logger.info('LS Instantiate and train LASSO: '
                    'CV alpha_weight hyperparameter optimization...')

        kf = KFold(n_splits=cv_folds, random_state=rand_seed, shuffle=shuffle) # Define the split
        alpha_optimal = L1_alpha[0]
        r2_optimal = -np.inf
        for alpha_ in L1_alpha:
            logger.info(('alpha_ = %f') % alpha_)

            LS['r2_train'][alpha_] = []
            LS['r2_cv'][alpha_] = []
            i_cv = 0
            for train_index, cv_index in kf.split(rates):
                LS[alpha_+i_cv] = linear_model.Lasso(alpha=alpha_, 
                                                     fit_intercept=False,
                                                     random_state=rand_seed, 
                                                     selection='random',
                                                     max_iter=max_iter_) 
                train_features, cv_features = np.array(features.loc[train_index,:]), np.array(features.loc[cv_index,:])
                train_rates, cv_rates = np.array(rates[rr+'z'+ts+ot][train_index]), np.array(rates[rr+'z'+ts+ot][cv_index])

                logger.info('LASSO1 fold %d training...' % i_cv)
                LS[alpha_+i_cv].fit(train_features, train_rates)
                LS['r2_train'][alpha_].append(LS[alpha_+i_cv].score(train_features, train_rates))
                LS['r2_cv'][alpha_].append(LS[alpha_+i_cv].score(cv_features, cv_rates))

                logger.info('LASSO1 fold %d train R^2: %f' % (i_cv, LS['r2_train'][alpha_][-1])) 
                logger.info('LASSO1 fold %d CV R^2: %f' % (i_cv, LS['r2_cv'][alpha_][-1]))
                i_cv += 1 #Move to next CV fold

            logger.info('LASSO1 average train r2 %f' % np.mean(LS['r2_train'][alpha_]))
            r2_cv = np.mean(LS['r2_cv'][alpha_])
            logger.info('LASSO1 average CV r2 %f' % r2_cv)

            logger.info('LS Instantiate and train LASSO1 given alpha_ to get features')
            LS[alpha_] = linear_model.Lasso(alpha=alpha_, 
                                            fit_intercept=False,
                                            random_state=rand_seed, 
                                            selection='random',
                                            max_iter=max_iter_) 
            LS[alpha_].fit(np.array(features), np.array(rates[rr+'z'+ts+ot]))
            logger.info('Get LASSO1 feature importances = abs(coef)')
            LS['f'][alpha_] = pd.DataFrame()
            LS['f'][alpha_]['features'] = features.columns
            LS['f'][alpha_]['coef'] = LS[alpha_].coef_
            LS['f'][alpha_]['importance'] = abs(LS[alpha_].coef_)
            LS['f'][alpha_].sort_values(by='importance',ascending=False, inplace=True)
            filename = 'LASSO1_'+ts+'_'+rr+'_feat'+str(args.feats)+'_'+str(alpha_)+'_nonzero.tsv'
            logger.info('Write nonzero features to file %s' % filename) 
            LS['f'][alpha_][LS['f'][alpha_]['importance']>0].to_csv(os.path.join(path,'select_features',filename), 
                                                                    sep='\t', index=False)          

            if r2_cv > r2_optimal:
                alpha_optimal = alpha_
                r2_optimal = r2_cv

        LS['r2_train'] = pd.DataFrame.from_dict(LS['r2_train']) 
        LS['r2_cv'] = pd.DataFrame.from_dict(LS['r2_cv'])
        filename = 'LASSO1_'+ts+'_'+rr+'_feat'+str(args.feats)+'_r2_'
        LS['r2_train'].to_csv(os.path.join(path,'select_features',filename+'train.tsv'), sep='\t', index=False)       
        LS['r2_cv'].to_csv(os.path.join(path,'select_features',filename+'cv.tsv'), sep='\t', index=False) 
        logger.info('Write r2 train and CV to files %s[train|cv].tsv' % filename)


        logger.info('LASSO1 optimized alpha_weight %f' % alpha_optimal)
        logger.info('LASSO1 optimized CV r2 %f' % r2_optimal)                    

    
if __name__ == '__main__':
    main()