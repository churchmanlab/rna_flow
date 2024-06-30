#!/usr/bin/env python
# coding: utf-8

# ## fit.py
# Author: Robert Ietswaart, 20200418
# License: BSD2.  
# Python v3.7.4 
# 
# For Brendan's project: fit time scales to his subcellular timelapse seq using observed new to total ratio (Lambda)

import numpy as np
import new_total_ratio as ntr

def calc_res(k_fit, model, times, lam_obs, k_para=None):
    """Calculate Residual
       for input model
       lam_obs : observed new to total ratio lambda
       k_fit : parameter to fit in model
       k_para: other parameters in model
       times : list with observed time points >0
       model : function describing mathematical model
       """
    if model in [ntr.lam_chr,
                 ntr.lam_nucpl,
                 ntr.lam_nuc,
                 ntr.lam_nuc_from_chr,
                 ntr.lam_cyto,
                 ntr.lam_cyto_from_chr,
                 ntr.lam_total,
                 ntr.lam_poly,
                 ntr.lam_poly_from_chr,
                 ntr.lam_total,
                 ntr.lam_total_one_step,
                 ntr.lam_total_from_chr,
                 ntr.TCconv]:
        lam_predict = model(k_fit, times, k_para)
    elif model in [np.poly1d]:
        poly = model(k_fit)
        lam_predict = poly.__call__(times)
    res = lam_predict - lam_obs
    return res

def calc_res_2reps(k_fit, model, times, lam_obs_rep1, lam_obs_rep2, k_para=None):
    """Calculate Residual
       for input model
       lam_obs : observed new to total ratio lambda
       k_fit : parameter to fit in model
       k_para: other parameters in model
       times : list with observed time points >0
       model : function describing mathematical model
       """
    if model in [ntr.lam_nuc,
                 ntr.lam_cyto,
                 ntr.lam_total,
                 ntr.lam_poly,
                 ntr.TCconv]:
        lam_predict = model(k_fit, times, k_para)
    elif model in [np.poly1d]:
        poly = model(k_fit)
        lam_predict = poly.__call__(times)
    res1 = lam_predict - lam_obs_rep1
    res2 = lam_predict - lam_obs_rep2
    return res1.append(res2, ignore_index=True)

def calc_res_3var(k, model1, model2, model3, times, lam_obs):
    """Calculate 2d Residuals sum of squares (RSS)
       for input model
       lam_obs : observed new to total ratio lambda
       k : list with parameters of model
       times : list with observed time points >0
       model : function describing mathematical model
       """
    lam_predict = [model1(k, times), 
                   model2(k, times),
                   model3(k, times)]
    res = lam_predict - lam_obs
    return res

def calc_res_matrix(params, model, observations, fixed_params):
    return np.concatenate((model(params, observations, fixed_params) - observations), axis=None)

def _aic(n,k,rss):
    """first equation on page e7 of McShane et al
       k : number of parameters in model"""
    aic = 2.0*k + n*np.log(rss/n) + 2.0*k*(k+1)/(n-k-1)
    return aic
 
def calc_aic(n, rss1, rss2, k1, k2):
    """Calculate Akaike Information Criterium (AIC)
       from residuals for 1 and 2 state models
       n : number of data points
       rss1 : residual sum of squares of model 1
       rss2 : as rss1 for model 2
       k1 : number of free parameters of model 1
       k2 : as k1 for model 2
       """
    aic1 = _aic(n, k1, rss1) 
    aic2 = _aic(n, k2, rss2)
    aic_min = min(aic1, aic2)
    sum_of_probs = np.exp((aic_min - aic1) / 2) + np.exp((aic_min - aic2) / 2)
    prob_aic1 = np.exp((aic_min - aic1) / 2) / sum_of_probs
    prob_aic2 = np.exp((aic_min - aic2) / 2) / sum_of_probs
    return [prob_aic1, prob_aic2]

def get_tc_types_for_gene_jit(gene, r, fr, GS_keys, TB):
    suffix = '_below'
    TC_TYPES_gene = []
    if r+fr+'top1000' in GS_keys:
        TC_TYPES_gene = np.asarray([1,0], dtype='int32')
        if gene.startswith('ENSMUSG') or gene.startswith('ENSG'):
            id_type = 'ENS_ID'
        else:
            id_type = 'Symbol'
        
        if (gene in TB[r+fr+'top1000'][id_type].values):
            TC_TYPES_gene = np.asarray([1], dtype='int32')
        elif (gene in TB[r+fr+'bottom500'][id_type].values) or \
             (gene in TB[r+fr+'bottom500'+suffix][id_type].values): 
            TC_TYPES_gene = np.asarray([0], dtype='int32')
    return TC_TYPES_gene