#!/usr/bin/env python
# coding: utf-8

# ## RNA decay fitting of models
# Author: Robert Ietswaart  
# Date: 20210712  
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4, pytorch v1.4  
# 
# For Brendan's project: perform timescale fitting based on Grand-Slam new to total ratio MAP values for individual genes for human K562.  
# Source: `Timescale_fit_20210603.ipynb`  
# GS
# TC background rates by Brendan  
# TC conversion rates estimated with BM model (`TCconversion_from_background_20210625_human.ipynb`).


import os
import re
import pandas as pd
import numpy as np
import math
import copy
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns
from scipy.optimize import least_squares
from scipy.stats import beta, pearsonr
import new_total_ratio as ntr
import fit

np.random.seed(12345)



path = os.path.join('/n','groups','churchman','ri23','bseq','GS20210713_human') #'GS20210625_human') #
outpath = os.path.join('/n','groups','churchman','ri23','bseq','MLE20210712_human')

fracs = ['chr', 'nucpl', 'nuc', 'cyto', 'poly', 'tot']
fracs_model = {'chr': ['chr_fit'], 
               'nucpl': ['nucpl_fit'], 
               'nuc': ['nuc_pred_from_chr_nucpl', 'nuc_fit'],
               'cyto': ['cyto_fit_from_chr_nucpl', 'cyto_fit_from_nuc'],
               'poly': ['poly_fit_from_chr_nucpl', 'poly_fit_from_nuc'],
               'tot': ['tot_pred_from_chr_nucpl', 'tot_pred_from_nuc', 'tot_fit']}
reps = ['T','U']
time_id = [str(i) for i in range(1,6)]
background_id = {r: '1' for r in reps}
time_mins = [0,15,30,60,120]
time_measured = pd.Series(time_mins[1:])
T_max = time_mins[-1]+30 #only used for plotting continuous curves
time_cont = pd.Series(range(0, T_max)) #only used for plotting continuous curves

TC_TYPES = ['top1000','bottom500']#BM model turnover method to estimate TC
      
GS = dict()
TB = dict()#TopBottom genes

for r in reps:
    for fr in fracs:
        for tc in TC_TYPES:
            filename = r + '_' + fr + '_noMT_' + tc + '.csv'
            if os.path.exists(os.path.join(path, filename)):
                GS[r+fr+tc]= pd.read_csv(os.path.join(path, filename) ,index_col=0)
                
                if tc == 'top1000':
                    SUFFICES = ['genes.turnover']
                if tc == 'bottom500':
                    SUFFICES = ['', '_below']
                    
                for suffix in SUFFICES:
                    filename = r + '1-' + fr + '_' + tc + suffix + '.MAPs.txt'
                    
                    if tc == 'top1000':
                        TB[r+fr+tc]= pd.read_csv(os.path.join(path, filename), 
                                                 sep='\t', header=None, 
                                                 names=['ENS_ID', 'Symbol', 'MAP'])
                    if tc == 'bottom500':
                        TB[r+fr+tc+suffix]= pd.read_csv(os.path.join(path, filename), 
                                                 sep='\t', header=None, 
                                                 names=['ENS_ID', 'Symbol', 'MAP'])


def get_tc_types_for_gene(gene, r, fr):
    TC_TYPES_gene = []
    if r+fr+'top1000' in GS.keys():
        TC_TYPES_gene = copy.deepcopy(TC_TYPES)
        if gene.startswith('ENSG'):
            id_type = 'ENS_ID'
        else:
            id_type = 'Symbol'
        
        if (gene in TB[r+fr+'top1000'][id_type].values):
            TC_TYPES_gene.remove('bottom500')
        elif (gene in TB[r+fr+'bottom500'][id_type].values) or (gene in TB[r+fr+'bottom500'+suffix][id_type].values):
            TC_TYPES_gene.remove('top1000')   
    return TC_TYPES_gene




# ### Fitting: initialize rate vector
# 
# references:  
# scipy.optimize.least_squares https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
# 
# tutorial robust nonlinear regression:
# https://scipy-cookbook.readthedocs.io/items/robust_regression.html


N_para = 8 
eps = 1e-6
k_bounds = (np.zeros(N_para)+eps,np.asarray([np.inf for i in range(N_para)]))

# k0**(-1) = \
Timescales = ['T_release', 
              'T_export_from_chr', 
              'T_export_from_nuc',
              'T_cytodeg_from_chr',
              'T_cytodeg_from_nuc',
              'T_translation_from_chr', 
              'T_translation_from_nuc',
              'T_totdeg']

k0 = abs(np.random.normal(1, 0.3, N_para))


def get_model_fit_parts(r, frm, tc):
    fit_idx = -1 
    model = None
    fixed_para = None
    if frm == 'chr_fit':
        fit_idx = 0
        model = ntr.lam_chr
    elif frm == 'nucpl_fit':
        fit_idx = 1
        model = ntr.lam_nucpl
        fixed_para = k_fit[tc+r][0]
    elif frm == 'nuc_fit':
        fit_idx = 2
        model = ntr.lam_nuc
    elif frm == 'cyto_fit_from_chr_nucpl':
        fit_idx = 3
        model = ntr.lam_cyto_from_chr
        fixed_para = [k_fit[tc+r][0], k_fit[tc+r][1], k0[5]]
    elif frm == 'cyto_fit_from_nuc':
        fit_idx = 4
        model = ntr.lam_cyto
        fixed_para = k_fit[tc+r][2]
    elif frm == 'poly_fit_from_chr_nucpl':
        fit_idx = 5
        model = ntr.lam_poly_from_chr 
        fixed_para = [k_fit[tc+r][0], k_fit[tc+r][1], 
                      k_fit[tc+r][3]]
    elif frm == 'poly_fit_from_nuc':
        fit_idx = 6
        model = ntr.lam_poly
        fixed_para = [k_fit[tc+r][2], k_fit[tc+r][4]]
    elif frm == 'tot_fit':
        fit_idx = 7
        model = ntr.lam_total_one_step

    if fit_idx >= 0:
        init_para = k0[fit_idx:(fit_idx+1)]
        bounds_para = [k_bounds[0][fit_idx], k_bounds[1][fit_idx]]
    else:
        init_para = None
        bounds_para = None

    return model, init_para, bounds_para, fixed_para, fit_idx


def get_model_pred_para(r, frm, tc):
    if frm == 'chr_fit':
        model = ntr.lam_chr
        pred_para = k_fit[tc+r][:1]
    elif frm == 'nucpl_fit':
        model = ntr.lam_nucpl
        pred_para = k_fit[tc+r][:2]
    elif frm == 'nuc_fit':
        model = ntr.lam_nuc
        pred_para = k_fit[tc+r][2:3]
    elif frm == 'cyto_fit_from_chr_nucpl':
        model = ntr.lam_cyto_from_chr
        pred_para = [k_fit[tc+r][0], k_fit[tc+r][1],
                     k_fit[tc+r][5], k_fit[tc+r][3]]
    elif frm == 'cyto_fit_from_nuc':
        model = ntr.lam_cyto
        pred_para = [k_fit[tc+r][2], k_fit[tc+r][6],
                     k_fit[tc+r][4]]
    elif frm == 'poly_fit_from_chr_nucpl':
        model = ntr.lam_poly_from_chr
        pred_para = [k_fit[tc+r][0], k_fit[tc+r][1],
                     k_fit[tc+r][5], k_fit[tc+r][3]]
    elif frm == 'poly_fit_from_nuc':
        model = ntr.lam_poly
        pred_para = [k_fit[tc+r][2], k_fit[tc+r][6],
                     k_fit[tc+r][4]]
    elif frm == 'tot_fit':
        model = ntr.lam_total_one_step
        pred_para = k_fit[tc+r][7:]
    elif frm == 'nuc_pred_from_chr_nucpl':
        model = ntr.lam_nuc_from_chr
        pred_para = k_fit[tc+r][:2]
    elif frm == 'tot_pred_from_chr_nucpl':
        model = ntr.lam_total_from_chr
        pred_para = [k_fit[tc+r][0],k_fit[tc+r][1],
                     k_fit[tc+r][5],k_fit[tc+r][3]]
    elif frm == 'tot_pred_from_nuc':
        model = ntr.lam_total
        pred_para = [k_fit[tc+r][2],k_fit[tc+r][6],k_fit[tc+r][4]]

    return model, pred_para


def get_chi2(ntr_meas,ntr_model):
    eps = 1e-16
    chi2 = 0
    for i in ntr_meas.index:
        chi2 = chi2 + (ntr_model[i]-ntr_meas[i])**2/(ntr_model[i]+eps)
    return chi2



estimate_type = 'MAP' 

#determine genes that have sufficient data to estimate timescales 
#NB: top1000 and bottom500 GS have identical gene ids
genes_w_rates = dict()
for r in reps:
    for fr in fracs:
        if r+fr+'top1000' in GS.keys():
            if fr == 'chr':
                genes_w_rates[r] = set(GS[r+fr+'top1000']['Gene'])
                print(r, fr, len(genes_w_rates[r]))
            elif fr in ['nucpl','cyto','poly','nuc','tot']:
                genes_w_rates[r] = genes_w_rates[r].union(GS[r+fr+'top1000']['Gene'])
                print(r, fr, len(genes_w_rates[r]))
genes_w_rates = genes_w_rates['T'].union(genes_w_rates['U'])
genes_w_rates = sorted(list(genes_w_rates))
print('total number of genes: union of T and U:', len(genes_w_rates))



#initialize the fit dictionary with columns
fits = dict()
fits['Gene'] = []
fits['Symbol'] = []
for r in reps:
    for tc in TC_TYPES:
        for timescale in Timescales:
            fits[r+'.'+tc+'.'+timescale] = []
for r in reps:
    for tc in TC_TYPES:
        for fr in fracs:
            if r+fr+tc in GS.keys():
                for frm in fracs_model[fr]:
                    fits[r+'.'+tc+'.'+frm+'.chi2'] = []

#Fit all genes
count = 0
for ensid in genes_w_rates:  
    if (count % 100) == 0:
        print(count,'/',len(genes_w_rates))
    
    #initialize k_fit per gene
    k_fit = dict()
    NTR = dict()
    NTR_model = dict() #NTR_pred_measured
    for r in reps:
        for tc in TC_TYPES:
            k_fit[tc+r] = np.asarray([np.nan for t in Timescales])
                
                
    for r in reps:
        gs_times = [r+t+'.'+estimate_type for t in time_id[1:]]
        for fr in fracs:
            #check which TC conversion rate is applicable to the gene
            TC_TYPES_gene = get_tc_types_for_gene(ensid, r, fr) 
#             print(TC_TYPES_gene)
            for tc in TC_TYPES_gene:
                if r+fr+tc in GS.keys():
                    try:
                        gene_idx = GS[r+fr+tc][GS[r+fr+tc]['Gene']==ensid].index[0]
                        symbol = GS[r+fr+tc]['Symbol'][gene_idx]
                        NTR[r+fr+tc] = pd.Series(list(GS[r+fr+tc].loc[gene_idx,gs_times])) 

                        for frm in fracs_model[fr]:
                            model, init_para, bounds_para, fixed_para, fit_idx = get_model_fit_parts(r, frm, tc)

                            if model is not None:
                                try:
                                    fit_val = least_squares(fit.calc_res, 
                                                            init_para, 
                                                            args=(model, time_measured, NTR[r+fr+tc], fixed_para),
                                                            bounds=bounds_para,
                                                            gtol=1e-14, ftol=1e-14,
                                                            loss='linear')

                                    k_fit[tc+r][fit_idx] = fit_val.x[0]

                                    if len(TC_TYPES_gene) == 1:#apply this fitted rate to both top and bottom
                                        if tc == 'top1000':
                                            tc2 = 'bottom500'
                                        elif tc == 'bottom500':
                                            tc2 = 'top1000'
                                        k_fit[tc2+r][fit_idx] = fit_val.x[0]
                                except (ValueError):
                                    pass
                    except (IndexError, KeyError):
                        pass


        #append rates in fits dictionary            
        if r in reps:
            for tc in TC_TYPES:
                for i, timescale in enumerate(Timescales):
                    fits[r+'.'+tc+'.'+timescale].append(1/k_fit[tc+r][i])
        
        #Get predictions
        for fr in fracs:           
            for tc in TC_TYPES:
                if r+fr+tc in GS.keys():
                    try:
                        gene_idx = GS[r+fr+tc][GS[r+fr+tc]['Gene']==ensid].index[0]
                        symbol = GS[r+fr+tc]['Symbol'][gene_idx]
                        NTR[r+fr+tc] = pd.Series(list(GS[r+fr+tc].loc[gene_idx,gs_times])) 

                        for frm in fracs_model[fr]:                 
                            model, pred_para = get_model_pred_para(r, frm, tc)
                            NTR_model[r+frm+tc] = model(pred_para, time_measured)

                            chi2 = get_chi2(NTR[r+fr+tc], NTR_model[r+frm+tc])
                            fits[r+'.'+tc+'.'+frm+'.chi2'].append(chi2)  
                    except (IndexError,KeyError):
                        for frm in fracs_model[fr]:
                            fits[r+'.'+tc+'.'+frm+'.chi2'].append(np.nan)
                    
    fits['Gene'].append(ensid)    
    fits['Symbol'].append(symbol)
            
    count+=1

fits = pd.DataFrame(fits)
filename = 'MLE20210804_human_fit_timescales.tsv'
fits.to_csv(os.path.join(outpath, filename), sep='\t',index=False)


