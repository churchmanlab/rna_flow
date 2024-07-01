#!/usr/bin/env python
# coding: utf-8

# ## Subcellular Timelapse seq: TC conversion rate estimation
# Author: Robert Ietswaart  
# Date: 20231125  
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4, pytorch v1.4 
# 
# source: `TCconversion_from_background_20210625_human.ipynb`  
# Estimate with binomial mixture model the TC conversion rates from the full no4sU and treated sample distributions of human (K562) experiments `T` and `U` for all available fractions.

import os
import pandas as pd
import numpy as np
import math
import copy
import logging
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns
from lmfit import minimize, Parameters

from __init__ import default_logger_format, default_date_format

import new_total_ratio as ntr
import fit

def main(): 

    outpath = os.path.join('/n','groups','churchman','ri23','bseq','TC20231125_K562')
    reps = ['T','U']
    background_id = {r: '1' for r in reps}
    time_id = [str(i) for i in range(1,6)]
    fracs = ['chr','nuc', 'cyto', 'poly','tot']
    time_mins = [0, 15, 30, 60, 120]
    time_measured = pd.Series(time_mins[1:])
    path = dict() 
    path['chr'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2021-05-16_T_U', 'STAR_2023')
    path['nuc'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2021-04-07_T_U', 'STAR_2023')
    path['cyto'] = path['chr']
    path['poly'] = path['chr']
    path['tot'] = path['nuc']
    TC_TYPES = ['top1000genes_turnover', 'bottom500genes_turnover']

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('TC')
    log_file = os.path.join(outpath,'LogErr', 'TCconversion_from_background_20231125_K562_py.log')
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)


    RAW = dict()
    A = dict()
    for r in reps:
        for fr in fracs:
            for tb in TC_TYPES:
                for t in time_id:
                    if tb == 'top1000genes_turnover':
                        path_data = os.path.join(path[fr], r + t + '_' + fr,
                                                 r + t + '_' + fr + '_MM_temp_turnover')
                    elif tb == 'bottom500genes_turnover':
                        path_data = os.path.join(path[fr], r + t + '_' + fr,
                                                 r + t + '_' + fr + '_MM_temp_turnover_slow')

                    filename = r + t + '_' + fr + '_'+ tb +'_TcountANDTCconv.txt'

                    if os.path.exists(os.path.join(path_data, filename)):
                        RAW[r+fr+tb+t]= pd.read_csv(os.path.join(path_data, filename), sep='\t', index_col=0)

                        #creating A[n,k] matrix: insert zeros where RAW has not missing lines/columns
                        N = RAW[r+fr+tb+t].index[-1] + 1
                        K = int(RAW[r+fr+tb+t].columns[-1]) + 1
                        A[r+fr+tb+t] = np.zeros((N,K))
                        for n in RAW[r+fr+tb+t].index:
                            for k in RAW[r+fr+tb+t].columns:
                                A[r+fr+tb+t][n, int(k)] = RAW[r+fr+tb+t][k][n]


    out = dict()
    for r in reps:
        for fr in fracs:
            for tb in TC_TYPES:
                if r+fr+tb+background_id[r] in A.keys():
                    #Fit 1 error rates background
                    params1 = Parameters()
                    params1.add('Error_1', value=0.0001, min=0, max=1)
                    out[r+fr+tb+'err1'] = minimize(fit.calc_res_matrix, params1, method='leastsquare',
                                             args=(ntr.one_error_background, A[r+fr+tb+background_id[r]], None))

                    #Fit 2 error rates background
                    params2 = Parameters()
                    params2.add('Error_1', value=0.0001, min=0, max=1)
                    params2.add('Error_2', value=0.00001, min=0, max=1)
                    params2.add('frac_err', value=0.5, min=0, max=1)
                    out[r+fr+tb+'err2'] = minimize(fit.calc_res_matrix, params2, method='leastsquare',
                                                args=(ntr.two_error_background, A[r+fr+tb+background_id[r]], None))

                    #Confirm with AIC that 2 error model better explains the data
                    prob1, prob2 = fit.calc_aic(out[r+fr+tb+'err1'].ndata,
                                                out[r+fr+tb+'err1'].chisqr,
                                                out[r+fr+tb+'err2'].chisqr,
                                                out[r+fr+tb+'err1'].nvarys,
                                                out[r+fr+tb+'err2'].nvarys)
                    logger.info('%s %s %s probability that model1 is favored %f' % (r, fr, tb, prob1))
                    logger.info('%s %s %s probability that model2 is favored %f' % (r, fr, tb, prob2))


    # ### Fit pc the TC conversion rate on treated samples
    TC_INITS = [0.001,0.01]
    FRAC_INITS = [0.05, 0.5]

    paramsTC = Parameters()
    for r in reps:
        for fr in fracs:
            logger.info('%s %s' % (r,fr))
            for tb in TC_TYPES:
                logger.info(tb)
                for t in time_id:
                    logger.info(t)
                    if (not (t == background_id[r])) and (r+fr+tb+t in A.keys()): 
                        chi2 = -1
                        for tc_init in TC_INITS:
                            for frac_init in FRAC_INITS:
                                paramsTC.add('TC', value=tc_init, min=0, max=1)
                                paramsTC.add('frac', value=frac_init, min=0, max=1)
                                out_temp = minimize(fit.calc_res_matrix, paramsTC, #method='leastsq',
                                            args=(ntr.TC_rate_from_background, 
                                                  A[r+fr+tb+t], 
                                                  out[r+fr+tb+'err2'].params))
                                if (chi2 < 0) or (chi2 > out_temp.chisqr):
                                    out[r+fr+tb+t] = out_temp
                                    chi2 = out_temp.chisqr


    # # Write the predicted rates to file
    for r in reps:
        for fr in fracs:
            for tb in TC_TYPES:
                if r+fr+tb+background_id[r] in A.keys():
                    output = pd.DataFrame()
                    output['time'] = time_measured
                    output_params = ['TC', 'frac']
                    for para in output_params:
                        output[para] = [out[r+fr+tb+t].params[para].value for t in time_id[1:]]
                        output[para+'_stderr'] = [out[r+fr+tb+t].params[para].stderr for t in time_id[1:]]
                    filename = r+'_'+fr+'_'+tb+'_BMmodel_TC_rate_predictions' 
                    output.to_csv(os.path.join(outpath, filename + '.tsv'), sep='\t', index=False)
                    #also write full output: for completeness
                    with open(os.path.join(outpath, filename + '_full.json'), 'w') as outfile:
                        out[r+fr+tb+t].params.dump(outfile)
                    logger.info('%s %s %s' % (r, fr, tb))
                    logger.info(output.head())

                    output = pd.DataFrame()
                    output['time'] = [0]
                    output_params = ['Error_1', 'Error_2', 'frac_err']
                    for para in output_params:
                        output[para] = [out[r+fr+tb+'err2'].params[para].value]
                        output[para+'_stderr'] = [out[r+fr+tb+'err2'].params[para].stderr]
                    filename = r+'_'+fr+'_'+tb+'_BMmodel_background_rate_predictions'
                    output.to_csv(os.path.join(outpath, filename + '.tsv'), sep='\t', index=False)
                    #also write full output: for completeness
                    with open(os.path.join(outpath, filename + '_full.json'), 'w') as outfile:
                        out[r+fr+tb+'err2'].params.dump(outfile)
                    logger.info(output.head())
                else: 
                    logger.info('no need to export rates for %s %s %s  since no input files' % (r, fr,tb))

if __name__ == '__main__':
    main()