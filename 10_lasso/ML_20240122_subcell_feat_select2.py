#!/usr/bin/env python
# coding: utf-8

# ## Machine learning to explain rates
# Author: Robert Ietswaart  
# Date: 20240122  
# License: BSD2. 
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4
# 
# Source: `ML_20230329_subcell_feat_select2.py`  
# For RNA flow project.



import os
import re
import copy
import numpy as np
import pandas as pd
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
import logging
import argparse
from scipy.stats import ttest_1samp, ttest_rel #ttest_ind,

from sklearn.model_selection import train_test_split, KFold
from sklearn import linear_model

from __init__ import default_logger_format, default_date_format 



def main():
    np.random.seed(12345)

    parser = argparse.ArgumentParser(
        description='Machine learning feature selection round 2 and test on Bayes rates.')
    parser.add_argument('--rate',type=str, default='whole_cell',help='provide a rate type. Choose from: '
                        'chr, chr_release, nuc, nucexp, cyto, poly_entry, whole_cell')

    num2feat_files = {1: 'gene_structure_v2.txt gene_structure_agarwal.txt',
                      2: 'gene_sequence_v2.txt gene_sequence_agarwal.txt',
                      3: 'Codon_features.tsv.gz', 4: 'kmer_ORF_features.tsv.gz',
                      5: 'kmer_5UTR_features.tsv.gz', 6: 'kmer_3UTR_features.tsv.gz',
                      7: 'polyA_tail_length_data.txt', 8: 'RNA_modifications_v2.txt',
                      9: 'gene_location_v2.txt', 10: 'histone_modifications_v1.txt',
                      11: 'microRNA_targets_v2.txt', 12: 'RBP_targets_filter_v2.txt',
                      13: 'TF_scRNAseq_v1.csv', 14: 'TF_MSigDB_targets_v2.tsv'}

    args = parser.parse_args()

    path = os.path.join('/n','groups','churchman','ri23','bseq','ML20240122')
    feature_path = os.path.join('/n','groups','churchman','ri23','bseq','RF20220426','features')

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('ML')
    log_file = os.path.join(path,'LogErr', 'ML_20240122_subcell_select2_k_%s_py.log' % args.rate)
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

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
    path_b = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240120_K562')
    filename_b = 'Bayes_Rates_20240120_'+ org_map[o] + '_final.tsv'
    path_k = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor20240112')
    filename_k = 'Bayes_factor_20240112_' + org_map[o] + '_final.tsv'

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
            C[o][rr+'.'+ts+ot] = copy.deepcopy(C[o][rr+'.'+ts+'_from_nucdeg.MAP'].where(
                    C[o]['PUND'], C[o][rr+'.'+rt+'chr'+ot]))
            ts = rt + 'nucdeg'#only for nucdeg genes according to Bayes Factor
            C[o][rr+'.'+ts+ot] = copy.deepcopy(C[o][rr+'.'+ts+'.MAP'].where(C[o]['PUND'], np.nan))
            ts = rt + 'nucexp'
            C[o][rr+'.'+ts+ot] = copy.deepcopy(C[o][rr+'.'+ts+'_from_nucdeg.MAP'].where(
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
    #identical split as in LASSO1 round
    C_cv = C[o].loc[train_idx,:]
    C_test = C[o].loc[test_idx,:]



    ts = rt + args.rate
    ot = OUT_TYPES[0]
    logger.info('Load LASSO round 1 CV results on rate %s' % ts+ot)

    RND1 = dict()

    for fk in num2feat_files.keys():
        fk = str(fk)
        for rr in org_red_reps[o]:
            filename = 'LASSO1_'+ts+'_'+rr+'_feat'+fk+'_r2_' 
            logger.info('rep %s feature set %s' % (rr,fk))
            RND1[fk+rr+'r2_train'] = pd.read_csv(os.path.join(path,'select_features',filename+'train.tsv'),sep='\t')
            RND1[fk+rr+'r2_cv'] = pd.read_csv(os.path.join(path,'select_features',filename+'cv.tsv'),sep='\t') 

    logger.info('t-tests to identify optimal alpha hyperparameter given CV results \n %s', ts)

    alpha_sig = 0.05

    params_opt = dict() #optimal parameters
    params_opt['feat'] = []
    params_opt['feat_files'] = []
    params_opt['alpha'] = []
    params_opt['pval'] = []
    params_opt['r2_cv'] = []
    for fk in num2feat_files.keys():

        fk = str(fk)
        params_opt['feat'].append(fk)
        params_opt['feat_files'].append(num2feat_files[int(fk)])
        L1_alpha = RND1[fk+org_red_reps[o][0]+'r2_cv'].columns.values

        pval_opt = 1
        alpha_opt = np.nan
        for a in L1_alpha:
            test_samples = []
            for rr in org_red_reps[o]:
                test_samples.extend(RND1[fk+rr+'r2_cv'][a].values)
            if pval_opt == 1: #model is not yet better than 0: do 1 sample test 
                stat, pval = ttest_1samp(test_samples, popmean=0,nan_policy='omit',alternative='greater')
                r2_cv_opt = np.mean(test_samples)
            else:#paired two-sample test
                stat, pval = ttest_rel(test_samples, test_samples_null,nan_policy='omit',alternative='greater')
            if pval < alpha_sig:#better than best hypothesis so far
                stat, pval = ttest_1samp(test_samples, popmean=0,nan_policy='omit',alternative='greater')
                pval_opt = pval
                alpha_opt = a 
                test_samples_null = test_samples
                r2_cv_opt = np.mean(test_samples)
        params_opt['alpha'].append(alpha_opt)
        params_opt['pval'].append(pval_opt)
        params_opt['r2_cv'].append(r2_cv_opt)

    params_opt = pd.DataFrame.from_dict(params_opt)
    logger.info(params_opt)
    r2_sum = np.sum(params_opt[~params_opt['alpha'].isna()]['r2_cv'])
    logger.info('sum of R2 (r) CV of individual relevant feature sets: %f (%f) \n' 
                '= upper limit on final model R2' % (r2_sum, np.sqrt(r2_sum)))



    logger.info('Get reproducible Lasso round 1 nonzero features for optimal alphas from files \n %s' % ts) 

    n_features = []
    for fk in num2feat_files.keys():
        fk = str(fk)
        logger.info('feature set %s' % fk)
        a_opt = params_opt[params_opt['feat']==fk]['alpha'].values[0]
        for rr in org_red_reps[o]:

            if type(a_opt) == str:
                filename = 'LASSO1_'+ts+'_'+rr+'_feat'+fk+'_'+str(a_opt)+'_nonzero.tsv'
                RND1[fk+rr+'f'] = pd.read_csv(os.path.join(path,'select_features',filename),sep='\t')
                if rr == org_red_reps[o][0]:
                    RND1[fk+'f'] = set(RND1[fk+rr+'f']['features'])
                else:
                    RND1[fk+'f'] = RND1[fk+'f'].intersection(set(RND1[fk+rr+'f']['features']))
            else:#a_opt = nan
                RND1[fk+rr+'f'] = set()
                RND1[fk+'f'] = set()
            logger.info('rep %s #features: %d' % (rr, len(RND1[fk+rr+'f'])))
        nf = len(RND1[fk+'f'])
        logger.info('reproducible #features: %d' % nf)
        n_features.append(nf)
        logger.info(RND1[fk+'f'])
    params_opt['n_features'] = n_features
    filename = 'LASSO2_'+ts+'_params_opt.tsv'
    params_opt.to_csv(os.path.join(path,'select_features',filename), sep='\t', index=False)        
    logger.info('Write params_opt to file %s' % filename)


    logger.info('Load relevant feature files %s:' % ts)

    F_raw = dict()

    feature_files = []
    for fk in num2feat_files.keys():
        a_opt = params_opt[params_opt['feat']==str(fk)]['alpha'].values[0]
        if type(a_opt) == str:#a_opt present, not nan
            feature_files.extend(num2feat_files[fk].split())       
    logger.info(feature_files)

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


    for fk in num2feat_files.keys():
        fk = str(fk)
        for f in num2feat_files[int(fk)].split():
            if f in feature_files:#only relevant features
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

                    feats = copy.deepcopy(F_raw[f][cols_merge])

                    #subset for relevant features
                    RND1[fk+'f'].add('Gene')
                    cols_rel = feats.columns[feats.columns.isin(RND1[fk+'f'])]
                    feats = feats[cols_rel]

                    #standardize quantitative features (get Z-score) to enhance learning        
                    feats[cols_rel[1:]] = (feats[cols_rel[1:]]-feats[cols_rel].mean()) / feats[cols_rel].std()

                elif f in {'gene_location_v2.txt', 'histone_modifications_v1.txt',
                           'microRNA_targets_v2.txt', 'RBP_targets_filter_v2.txt',
                           'TF_scRNAseq_v1.csv', 'TF_MSigDB_targets_v2.tsv'}:
                    logger.info('%s pivot, one hot encode %s' % (f, cols_merge))#no standardization needed

                    feats = F_raw[f][cols_merge].pivot(index=cols_merge[0], 
                                                       columns=cols_merge[1], 
                                                       values=cols_merge[1:])

                    for m, c in enumerate(cols_merge[1:]):
                        if m == 0:
                            feats[c] = feats[c].where(feats[c].isna(), 1)
                        feats[c] = feats[c].where(~feats[c].isna(), 0)
                    feats.columns = ['_'.join(col) for col in feats.columns.values] #from multi-index to index     
                    feats.reset_index(inplace=True) #Gene no longer index

                    feats = C[o][[sm_key]].merge(feats, on=sm_key, how ='left', suffixes=('', f))#feat row allrates

                    #subset for relevant features
                    RND1[fk+'f'].add('Gene')
                    cols_rel = feats.columns[feats.columns.isin(RND1[fk+'f'])]
                    feats = feats[cols_rel]
                else:
                    logger.info('fail: %s not recognized' % f)

                if f == feature_files[0]:#first file   
                    F = copy.deepcopy(feats)
                    logger.info('%s initial feature mat: #features: %d ,#genes: %d' % 
                                (f, len(F.columns)-1, len(F))) 
                else:
                    F = F.merge(feats, on=sm_key, how ='outer', suffixes=('', f))
                    logger.info('%s past merge feature mat: #features: %d ,#genes: %d' % 
                                (f, len(F.columns)-1, len(F)))


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

    logger.info('Merge rates and features of both bio reps')

    r1 = org_red_reps[o][0]
    r2 = org_red_reps[o][1]
    genes_w_rates1 = ~C_cv[r1+'z'+ts+ot].isna()
    genes_w_rates2 = ~C_cv[r2+'z'+ts+ot].isna()

    rates = pd.concat([ C_cv[genes_w_rates1][['Gene', r1+'z'+ts+ot]].rename(mapper={r1+'z'+ts+ot:ts}, axis=1),
                        C_cv[genes_w_rates2][['Gene', r2+'z'+ts+ot]].rename(mapper={r2+'z'+ts+ot:ts}, axis=1) ],
                      ignore_index=True)

    required_features = {'gene_structure_v2.txt', 'gene_structure_agarwal.txt', 'gene_sequence_v2.txt',
                        'gene_sequence_agarwal.txt', 'Codon_features.tsv.gz', 'kmer_ORF_features.tsv.gz',
                        'kmer_5UTR_features.tsv.gz', 'kmer_3UTR_features.tsv.gz', 'polyA_tail_length_data.txt', 
                        'RNA_modifications_v2.txt'}
    for f in feature_files:
        if f in required_features:
            rates = rates[rates['Gene'].isin(F_raw[f].dropna()['Gene'])]
    rates.sort_values(by=['Gene'], inplace=True)
    rates.reset_index(inplace=True, drop=True)

    features = pd.concat([ F[F['Gene'].isin(C_cv[genes_w_rates1]['Gene'])],
                           F[F['Gene'].isin(C_cv[genes_w_rates2]['Gene'])] ],
                         ignore_index=True)
    for f in feature_files:
        if f in required_features:
            features = features[features['Gene'].isin(F_raw[f].dropna()['Gene'])]
    features.sort_values(by=['Gene'], inplace=True)
    features.reset_index(inplace=True, drop=True)

    check1 = sum(rates['Gene'] == features['Gene'])
    logger.info('Confirm that rates and features are same order: '
                '%d = %d' % (len(rates), check1))
    rates.drop(['Gene'], axis=1, inplace=True)
    features.drop(['Gene'], axis=1, inplace=True) 

    L1_alpha = [0.0001, 0.00033, 0.00066, 0.001, 0.0033, 0.0066, 0.01, 0.033, 0.066, 0.1]
    
    n_cores = 1
    cv_folds = 10
    max_iter_ = 5e4 #default 1e3 (take 3mins), 1e6 can take >1.5h 
    logger.info('n LASSO round 2 CV folds: %d' % cv_folds)

    logger.info('LS Instantiate and train LASSO2: '
                'CV alpha_weight hyperparameter optimization...')

    LS = dict()
    #Lasso for feature selection
    LS['r2_train'] = dict()
    LS['r2_cv'] = dict()
    LS['r2_test'] = dict()
    LS['f'] = dict()

    kf = KFold(n_splits=cv_folds, random_state=rand_seed, shuffle=shuffle) # Define the CV split
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
            train_features, cv_features = features.loc[train_index,:], features.loc[cv_index,:]
            train_rates, cv_rates = rates[ts][train_index], rates[ts][cv_index]


            logger.info('LASSO2 CV fold %d training...' % i_cv)
            LS[alpha_+i_cv].fit(train_features, train_rates)
            LS['r2_train'][alpha_].append(LS[alpha_+i_cv].score(train_features, train_rates))
            LS['r2_cv'][alpha_].append(LS[alpha_+i_cv].score(cv_features, cv_rates))

            logger.info('LASSO2 fold %d train R^2: %f' % (i_cv, LS['r2_train'][alpha_][-1])) 
            logger.info('LASSO2 fold %d CV R^2: %f' % (i_cv, LS['r2_cv'][alpha_][-1]))
            i_cv += 1 #Move to next CV fold

        logger.info('LASSO2 average train r2 %f' % np.mean(LS['r2_train'][alpha_]))
        r2_cv = np.mean(LS['r2_cv'][alpha_])
        logger.info('LASSO2 average CV r2 %f' % r2_cv)

        logger.info('LS Instantiate and train LASSO2 given alpha_ on all training data to get features')
        LS[alpha_] = linear_model.Lasso(alpha=alpha_, 
                                        fit_intercept=False,
                                        random_state=rand_seed, 
                                        selection='random',
                                        max_iter=max_iter_) 
        LS[alpha_].fit(features, rates[ts])
        logger.info('Get LASSO2 feature importances = abs(coef)')
        LS['f'][alpha_] = pd.DataFrame()
        LS['f'][alpha_]['features'] = features.columns
        LS['f'][alpha_]['coef'] = LS[alpha_].coef_
        LS['f'][alpha_]['importance'] = abs(LS[alpha_].coef_)
        LS['f'][alpha_].sort_values(by='importance',ascending=False, inplace=True)
        filename = 'LASSO2_'+ts+'_alpha'+str(alpha_)+'_nonzero_feats.tsv'
        logger.info('Write nonzero features to file %s' % filename) 
        LS['f'][alpha_][LS['f'][alpha_]['importance']>0].to_csv(os.path.join(path,'select_features',filename), 
                                                                sep='\t', index=False)          

        if r2_cv > r2_optimal:
            alpha_optimal = alpha_
            r2_optimal = r2_cv

    LS['r2_train'] = pd.DataFrame.from_dict(LS['r2_train']) 
    LS['r2_cv'] = pd.DataFrame.from_dict(LS['r2_cv'])
    filename = 'LASSO2_'+ts+'_r2_'
    LS['r2_train'].to_csv(os.path.join(path,'select_features',filename+'train.tsv'), sep='\t', index=False)       
    LS['r2_cv'].to_csv(os.path.join(path,'select_features',filename+'cv.tsv'), sep='\t', index=False) 
    logger.info('Write r2 train and CV to files %s[train|cv].tsv' % filename)


    logger.info('LASSO2 optimized alpha_weight %f' % alpha_optimal)
    logger.info('LASSO2 optimized CV r2 %f' % r2_optimal)                    

    
    logger.info('Test performance: on unseen rates, merge of both bio reps')
    LS['r2_test'] = dict()

    genes_w_rates1 = ~C_test[r1+'z'+ts+ot].isna()
    genes_w_rates2 = ~C_test[r2+'z'+ts+ot].isna()

    rates = pd.concat([ C_test[genes_w_rates1][['Gene', r1+'z'+ts+ot]].rename(mapper={r1+'z'+ts+ot:ts}, axis=1),
                        C_test[genes_w_rates2][['Gene', r2+'z'+ts+ot]].rename(mapper={r2+'z'+ts+ot:ts}, axis=1) ],
                      ignore_index=True)
    for f in feature_files:
        if f in required_features:
            rates = rates[rates['Gene'].isin(F_raw[f].dropna()['Gene'])]
    rates.sort_values(by=['Gene'], inplace=True)
    rates.reset_index(inplace=True, drop=True)

    features = pd.concat([ F[F['Gene'].isin(C_test[genes_w_rates1]['Gene'])],
                           F[F['Gene'].isin(C_test[genes_w_rates2]['Gene'])] ],
                         ignore_index=True)
    for f in feature_files:
        if f in required_features:
            features = features[features['Gene'].isin(F_raw[f].dropna()['Gene'])]
    features.sort_values(by=['Gene'], inplace=True)
    features.reset_index(inplace=True, drop=True)

    check1 = sum(rates['Gene'] == features['Gene'])
    logger.info('Confirm that rates and features are same order: '
                '%d = %d' % (len(rates), check1))
    rates.drop(['Gene'], axis=1, inplace=True)
    features.drop(['Gene'], axis=1, inplace=True) 

    for alpha_ in L1_alpha:
        r2_test = LS[alpha_].score(features, rates[ts])
        logger.info('LASSO2  alpha %f on unseen test set r2 %f' % (alpha_, r2_test))

        LS['r2_test'][alpha_] = [r2_test]
    LS['r2_test'] = pd.DataFrame.from_dict(LS['r2_test']) 
    filename = 'LASSO2_'+ts+'_r2_test_seed'+str(rand_seed)+'.tsv'
    LS['r2_test'].to_csv(os.path.join(path,'select_features',filename), sep='\t', index=False)        
    logger.info('Write r2 test to file %s' % filename)
    logger.info('End')

    if args.rate == 'whole_cell':
        filename = 'Biorxiv2022Agarwal_half_lives.csv'
        agar_path = os.path.join('/n','groups','churchman','ri23','bseq','RF20220426')
        agar = pd.read_csv(os.path.join(agar_path, filename))

        agar.drop(['Symbol'], axis=1, inplace=True)
        agar.head()
        ts = RATE_TYPE[0] + 'whole_cell' #half lives
        
        logger.info('Agarwal: Master Train/CV vs %s test set split' % test_size)

        agar_cv, agar_test = train_test_split(agar, 
                                              test_size = test_size, 
                                              random_state = rand_seed, 
                                              shuffle = shuffle)
        
        logger.info('Agarwal rates, but out train/CV/test split')

        rates = copy.deepcopy(agar_cv)

        for f in feature_files:
            if f in required_features:
                rates = rates[rates['Gene'].isin(F_raw[f].dropna()['Gene'])]
        rates.sort_values(by=['Gene'], inplace=True)
        rates.reset_index(inplace=True, drop=True)

        features = copy.deepcopy(F[F['Gene'].isin(rates['Gene'])])
        for f in feature_files:
            if f in required_features:
                features = features[features['Gene'].isin(F_raw[f].dropna()['Gene'])]
        features.sort_values(by=['Gene'], inplace=True)
        features.reset_index(inplace=True, drop=True)

        check1 = sum(rates['Gene'] == features['Gene'])
        logger.info('Confirm that rates and features are same order: '
                    '%d = %d' % (len(rates), check1))
        rates.drop(['Gene'], axis=1, inplace=True)
        features.drop(['Gene'], axis=1, inplace=True) 
        
        logger.info('n Agarwal (rates only), our t/cv/t split folds: %d' % cv_folds)

        logger.info('LS Instantiate and train LASSO2 on Agarwal (rates only) with our train/CV/test split: '
                    'CV alpha_weight hyperparameter optimization...')

        LS = dict()
        #Lasso for feature selection
        LS['r2_train'] = dict()
        LS['r2_cv'] = dict()
        LS['r2_test'] = dict()
        LS['f'] = dict()

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
                train_features, cv_features = features.loc[train_index,:], features.loc[cv_index,:]
                train_rates, cv_rates = rates[ts][train_index], rates[ts][cv_index]


                logger.info('Agarwal (rates only), our t/cv/t split fold %d training...' % i_cv)
                LS[alpha_+i_cv].fit(train_features, train_rates)
                LS['r2_train'][alpha_].append(LS[alpha_+i_cv].score(train_features, train_rates))
                LS['r2_cv'][alpha_].append(LS[alpha_+i_cv].score(cv_features, cv_rates))

                logger.info('Agarwal (rates only), our t/cv/t split fold %d train R^2: %f' % (i_cv, LS['r2_train'][alpha_][-1])) 
                logger.info('Agarwal (rates only), our t/cv/t split fold %d CV R^2: %f' % (i_cv, LS['r2_cv'][alpha_][-1]))
                i_cv += 1 #Move to next CV fold

            logger.info('Agarwal (rates only), our t/cv/t split average train r2 %f' % np.mean(LS['r2_train'][alpha_]))
            r2_cv = np.mean(LS['r2_cv'][alpha_])
            logger.info('Agarwal (rates only), our t/cv/t split average CV r2 %f' % r2_cv)

            logger.info('LS Instantiate and train LASSO2 given alpha_ on all training data to get features')
            LS[alpha_] = linear_model.Lasso(alpha=alpha_, 
                                            fit_intercept=False,
                                            random_state=rand_seed, 
                                            selection='random',
                                            max_iter=max_iter_) 
            LS[alpha_].fit(features, rates[ts])
            logger.info('Get Agarwal (rates only), our t/cv/t split, feature importances = abs(coef)')
            LS['f'][alpha_] = pd.DataFrame()
            LS['f'][alpha_]['features'] = features.columns
            LS['f'][alpha_]['coef'] = LS[alpha_].coef_
            LS['f'][alpha_]['importance'] = abs(LS[alpha_].coef_)
            LS['f'][alpha_].sort_values(by='importance',ascending=False, inplace=True)
            filename = 'LASSO2_Agarwal_tcvtsplit_alpha'+str(alpha_)+'_nonzero_feats.tsv'
            logger.info('Write nonzero features to file %s' % filename) 
            LS['f'][alpha_][LS['f'][alpha_]['importance']>0].to_csv(os.path.join(path,'select_features',filename), 
                                                                    sep='\t', index=False)          

            if r2_cv > r2_optimal:
                alpha_optimal = alpha_
                r2_optimal = r2_cv

        LS['r2_train'] = pd.DataFrame.from_dict(LS['r2_train']) 
        LS['r2_cv'] = pd.DataFrame.from_dict(LS['r2_cv'])
        filename = 'LASSO2_Agarwal_tcvtsplit_r2_'
        LS['r2_train'].to_csv(os.path.join(path,'select_features',filename+'train.tsv'), sep='\t', index=False)       
        LS['r2_cv'].to_csv(os.path.join(path,'select_features',filename+'cv.tsv'), sep='\t', index=False) 
        logger.info('Write r2 train and CV to files %s[train|cv].tsv' % filename)
        logger.info('Agarwal rates: our train/CV/test split optimized alpha_weight %f' % alpha_optimal)
        logger.info('Agarwal rates: our train/CV/test split optimized CV r2 %f' % r2_optimal)  


        logger.info('Agarwal Test performance: on unseen rates')
        rates = copy.deepcopy(agar_test)

        for f in feature_files:
            if f in required_features:
                rates = rates[rates['Gene'].isin(F_raw[f].dropna()['Gene'])]
        rates.sort_values(by=['Gene'], inplace=True)
        rates.reset_index(inplace=True, drop=True)

        features = copy.deepcopy(F[F['Gene'].isin(rates['Gene'])])
        for f in feature_files:
            if f in required_features:
                features = features[features['Gene'].isin(F_raw[f].dropna()['Gene'])]
        features.sort_values(by=['Gene'], inplace=True)
        features.reset_index(inplace=True, drop=True)

        check1 = sum(rates['Gene'] == features['Gene'])
        logger.info('Confirm that rates and features are same order: '
                    '%d = %d' % (len(rates), check1))
        rates.drop(['Gene'], axis=1, inplace=True)
        features.drop(['Gene'], axis=1, inplace=True) 

        for alpha_ in L1_alpha:
            r2_test = LS[alpha_].score(features, rates[ts])
            logger.info('LASSO2  alpha %f on unseen test set r2 %f' % (alpha_, r2_test))

            LS['r2_test'][alpha_] = [r2_test]
        LS['r2_test'] = pd.DataFrame.from_dict(LS['r2_test']) 
        filename = 'LASSO2_Agarwal_tcvtsplit_r2_test_seed'+str(rand_seed)+'.tsv'
        LS['r2_test'].to_csv(os.path.join(path,'select_features',filename), sep='\t', index=False)        
        logger.info('Write r2 test to file %s' % filename)

    
if __name__ == '__main__':
    main()

