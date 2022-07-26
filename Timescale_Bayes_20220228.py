#!/usr/bin/env python
# coding: utf-8

# ## Bayesian inference to obtain rate posteriors
# Author: Robert Ietswaart  
# Date: 20220228  
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4

# Source: `Timescale_Bayes_20220222.ipynb`  
# For Subcellular Timelapse seq project: Bayesian inference to get timescale posteriors based on Grand-Slam new to total ratio posteriors for individual genes, now including **nuclear degradation**. 

# This script is used to run for all genes.  
# See also accompanying shell script `Timescale_Bayes_20220228.sh`.


import os
import numpy as np
import pandas as pd
import logging
import argparse
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns
from scipy.integrate import quad, dblquad, cumulative_trapezoid

# from __init__ import __version__
from __init__ import default_logger_format, default_date_format

import new_total_ratio_jit as ntr
import posteriors_jit as p
import fit

import numba as nb
from numba import jit
from numba.core import types
from numba.typed import Dict


def main():   
    np.random.seed(12345)
    batch_size = 100

    N_grid = 1000 #to divide [k_bound_lo,k_bound_hi] interval
    k_bound_lo = 1e-4 #unit: min^-1: 1 per 7 days
    k_bound_hi = 1e4 #unit: min^-1: 1 per 6 ms, if too restrictive: increase to 1e6
    eps_int = 1e-3 #numerical integration precision: default: 1e-3, can go low to 1e-6

    parser = argparse.ArgumentParser(
        description='Bayesian rate fitting on subcellular timelapse seq for standard and nucdeg model'
                    'in batches of %s genes.' % str(batch_size))
    parser.add_argument('--start_id',type=int, default=0,
                        help='row index of starting gene for batch.')

    args = parser.parse_args()

    path = os.path.join('/n','groups','churchman','ri23','bseq','GS20210506')
    outpath = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20220222')

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('Bayes_rate')
    log_file = os.path.join(outpath,'LogErr', 'scTLseq_bayes_20220228_batch%s_py.log' % str(args.start_id))
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)

    fracs = ['chr', 'nuc', 'cyto', 'poly', 'tot']
    fracs_model = {'chr': ['chr_fit'],  
                   'nuc': ['nuc_fit', 'nuc_fit_from_chr', 'nuc_fit_from_nucdeg'],
                   'cyto': ['cyto_fit_from_nucres', 'cyto_fit_from_chr', 'cyto_fit_from_nucdeg'],
                   'poly': ['poly_fit_from_nucres'],
                   'tot': ['tot_fit', 'tot_pred_from_nucres', 'tot_pred_from_chr', 
                           'tot_fit_from_nucdeg', 'tot_pred_from_nucdeg']}
    reps = ['R', 'S', 'G', 'H']
    red_r = {'G': 'G_R', 'H': 'H_S', 'R': 'G_R', 'S': 'H_S'} #reduced replicate to ensure full time series
    red_reps = ['G_R', 'H_S']#, 'T', 'U']
    red2reps = {'G_R': ['R','G'], 'H_S': ['S','H']}
    time_id = [str(i) for i in range(1,6)]
    background_id = {r: '1' for r in reps}
    time_mins = [0,15,30,60,120]
    time_measured = np.asarray(time_mins[1:])
    T_max = time_mins[-1]+30 #only used for plotting continuous curves
    time_cont = np.asarray(range(0, T_max)) #only used for plotting continuous curves


    TC_TYPES = ['top1000','bottom500']#BM model turnover method to estimate TC
    TC_from_jit = {1: 'top1000', 0: 'bottom500'}
    CI_para = ['alpha','beta']
    OUT_TYPES = ['.Mean', '.MAP', '.0.975.quantile', '.0.025.quantile'] # 
    #the CIs must match alpha/2 and 1-alpha/2 in get_post_ci()
    RATE_TYPE = ['half_life_','k_','T_']

    Timescales = ['T_chr',
                  'T_nuc',
                  'T_nucexp_from_chr',
                  'T_cyto',
                  'T_poly_entry',
                  'T_whole_cell',
                  'T_nucexp_from_nucdeg',
                  'T_nucdeg']


    Timescales2fracs = {'T_chr': ['chr'],
                        'T_nuc': ['nuc'],
                        'T_nucexp_from_chr': ['chr','nuc'],
                        'T_nucexp_from_nucdeg': ['tot'],
                        'T_nucdeg' : ['nuc'],
                        'T_cyto' : ['cyto'],
                        'T_poly_entry' : ['poly'],
                        'T_whole_cell' : ['tot']}

    GS = dict()
    TB = dict()#TopBottom genes
    for r in reps:
        for fr in fracs:
            for tc in TC_TYPES:
                filename = r + '_' + fr + '_noMT_' + tc + '.csv'
                if os.path.exists(os.path.join(path, filename)):
                    GS[r+fr+tc]= pd.read_csv(os.path.join(path, filename) ,index_col=0)

                    filename = r + '1-' + fr + '_' + tc + '.MAPs.txt'   
                    TB[r+fr+tc]= pd.read_csv(os.path.join(path, filename), 
                                             sep='\t', header=None, 
                                             names=['ENS_ID', 'Symbol', 'MAP'])
                    if tc == 'bottom500':
                        suffix = '_below'
                        filename = r + '1-' + fr + '_' + tc + suffix + '.MAPs.txt'
                        TB[r+fr+tc+suffix]= pd.read_csv(os.path.join(path, filename), 
                                                 sep='\t', header=None, 
                                                 names=['ENS_ID', 'Symbol', 'MAP'])


    # ### Posterior calculations

    #NB: NTR needs to be numba.typed.Dict and TC_TYPES_gene generated by calling fit.get_tc_types_for_gene_jit()
    @jit(nb.float64(nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), nopython=True) 
    def post_chr(kR, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_chr(kR, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_chr(kR, t)
        return prob


    @jit(nb.float64(nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), nopython=True) 
    def post_nuc(a, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_nuc(a, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_nuc(a, t)
        return prob        


    def post_nucexp_from_chr(kE, TC, NTR):
        '''
        TC and NTR are now dictionaries with chr and nuc fractions as keys.
        See N8p119 for derivation.
        '''
        def integrand_nucexp(z):
            return (post_nuc((kE * z * (kE + z)**(-1)), TC['nuc'], NTR['nuc']) * post_chr(z, TC['chr'], NTR['chr']))
        #integrate over small subdomains
        prob = 0
        prob_err = 0
        for i in range(8):
            k_int_lo = k_bound_lo * 10**(i) 
            k_int_hi = k_bound_lo * 10**(i+1)
            try:
                k_int, k_int_err = quad(integrand_nucexp, k_int_lo, k_int_hi, epsabs=eps_int)
    #             logger.info('post_nucexp_from_chr %d %f %f' % (i, k_int, k_int_err))
                if np.isnan(k_int):
                    logger.warning('post_nucexp_from_chr_err1: %f %f' % (k_int, k_int_err))
                    prob = np.nan
                    prob_err = np.nan
                    break
            except (IndexError, KeyError, ZeroDivisionError):
                logger.warning('post_nucexp_from_chr_err2')
                prob = np.nan
                prob_err = np.nan
                break
            prob += k_int
            prob_err += k_int_err
        return prob


    def post_nucdeg(ci_kE):
        '''
        Calculate the kN posterior by integration over kE with 
        ci_kE as integration boundaries.
        post: posterior over kN
        '''
        def marg_dist(kN_, TC_, NTR_):
            return quad(post_nuc, (kN_ + ci_kE[0]), (kN_ + ci_kE[1]), 
                        args=(TC_, NTR_), epsabs=eps_int)[0]
        post = marg_dist
        return post


    @jit(nb.float64(nb.float64, nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), 
         nopython=True) 
    def post_cyto(a, kD, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_cyto(a, kD, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_cyto_from_nuc(a, kD, t)
        return prob  


    @jit(nb.float64(nb.float64, nb.float64, nb.float64,
                    nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), nopython=True) 
    def post_poly_from_nuc(kNR, kD, kL, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_poly(kNR, kD, kL, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_poly_from_nuc(kNR, kD, kL, t)
        return prob 


    @jit(nb.float64(nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), nopython=True) 
    def post_tot(kDtot, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_total_one_step(kDtot, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_total_one_step(kDtot, t)
        return prob 


    @jit(nb.float64(nb.float64, nb.float64, nb.float64, 
                    nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), nopython=True) 
    def post_total_from_nucdeg(a, kD, kE, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_total_from_nucdeg(a, kD, kE, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_total_from_nucdeg(a, kD, kE, t)
        return prob


    def get_conditional_dist(post, TC, NTR, ci, nvar):
        '''
        Calculate the marginal posterior by integration over the previous rate with 
        ci as integration boundaries.
        post: multivariate posterior
        '''
        def marg_dist(x, TC, NTR):
            if nvar == 2:
                return quad(post, ci[0][0], ci[0][1], args=(x, TC, NTR), epsabs=eps_int)[0]
            elif nvar == 3:
                return dblquad(post, ci[0][0], ci[0][1], ci[1][0], ci[1][1], args=(x, TC, NTR), epsabs=eps_int)[0]
        cond_dist = marg_dist
        return cond_dist


    def get_post_mean(post, TC, NTR, norm_const):
        def integrand(x):#cannot be jitted bc it does not recognize function post
            return x * post(x, TC, NTR) / norm_const

        #integrate over small subdomains
        k_mean = 0
        k_mean_err = 0
        for i in range(8):#if k_bound_lo and k_bound_hi: update 4 or 8
            k_int_lo = k_bound_lo * 10**(i)
            k_int_hi = k_bound_lo * 10**(i+1)
            try:
                k_int, k_int_err = quad(integrand, k_int_lo, k_int_hi, epsabs=eps_int)
                logger.info('mean %f %f' % (k_int, k_int_err))
                if np.isnan(k_int):
                    logger.warning('mean_err1 %f %f' % (k_int_lo, k_int_hi))
                    k_mean = np.nan
                    k_mean_err = np.nan
                    break
            except (IndexError, KeyError, ZeroDivisionError):
                logger.warning('mean_err2 %f %f' % (k_int_lo, k_int_hi))
                k_mean = np.nan
                k_mean_err = np.nan
                break
            k_mean += k_int
            k_mean_err += k_int_err
        if k_mean == 0:
            logger.warning('mean_err3: mean is zero, set to k_bound_lo')
            k_mean = k_bound_lo#np.nan #eps_min #
            k_mean_err = k_bound_lo#np.nan #eps_min #       
        return k_mean, k_mean_err


    def get_post_ci(k_domain, post_grid):
        alpha = 0.05
        #CDF via cumulative trapezoid
        cdf = cumulative_trapezoid(post_grid, k_domain, initial=0) 
        omega = max(cdf)
        cdf = cdf / omega
        logger.info('norm_constant %f' % omega)
        k_temp = k_domain[cdf <= (alpha / 2)]
        if k_temp.size:
            ci_lo = max(k_temp)
        else:
            ci_lo = k_domain[0]
        k_temp = k_domain[cdf >= (1 - alpha / 2)]
        if k_temp.size:
            ci_hi = min(k_temp)
        else: 
            N_domain = len(k_domain)
            ci_hi = k_domain[N_domain -1]
        return ci_lo, ci_hi, omega

    
    def get_post_grid(post_, k_, TC_, NTR_):
        pg = np.empty(N_grid)
        for i in range(N_grid):
            pg[i] = post_(k_[i], TC_, NTR_) #evaluate posterior in a grid to get k_map and CIs
        mp = max(pg)
        return pg, mp  

    
    def increase_var_NTR(TC_, NTR_, vs):
        '''Generate a new posterior by adjusting NTR alphas and betas such that 
           the variance of beta_pdf is increased by factor vs but mean of beta_pdf 
           remains unaltered: math derivation on 20220227 notes'''
        if type(NTR_)==dict:
            for fr in NTR_.keys():
                NTR_[fr] = p.ntr_increased_var(TC_[fr], NTR_[fr], vs)
        else: 
            NTR_ = p.ntr_increased_var(TC_, NTR_, vs)
        return NTR_

    
    def get_estimates_from_post(post, TC, NTR, N_var, ci_k_prev=None):
        '''
        post : possibly multivariate posterior (function). If N_var > 1, the marginal
               distribution of the unknown rate, conditioned on ci_k_prev, will replace this
               multivariate posterior.
        N_var : number of rate arguments of the posterior.
        ci_k_prev : credible intervals of CIs of the previously fitted rates (list of tuples).
        ''' 
        if N_var > 1:
            post_new = get_conditional_dist(post, TC, NTR, ci_k_prev, N_var)
            post = post_new
        logger.info('g_e_f_p0')

        #Initialize the prior domain of rate
        k_domain = np.geomspace(k_bound_lo, k_bound_hi, num=N_grid)
        post_grid, max_post = get_post_grid(post, k_domain, TC, NTR)
        logger.info('max post %f' % max_post)
        logger.info('g_e_f_p1')

        if max_post == 0:
            logger.info('increase variance to get nonzero posterior') 
            for var_scale in range(2, 5):
                NTR_scale = increase_var_NTR(TC, NTR, var_scale)
                post_scale, max_scale = get_post_grid(post, k_domain, TC, NTR_scale)
                if max_scale > 0:
                    logger.info('max post increased: %f, var_scale %dx: %s' % (max_scale, var_scale, NTR_scale))
                    NTR = NTR_scale
                    post_grid = post_scale
                    max_post = max_scale
                    break
                else: 
                    logger.info('failed attempt to increase var:  var_scale %d' % var_scale)          

        if max_post > 0:
            map_k = k_domain[np.argmax(post_grid)]
            logger.info('g_e_f_p2') 

            ci_lo, ci_hi, norm_const = get_post_ci(k_domain, post_grid)
            post_grid = post_grid / norm_const
            logger.info('g_e_f_p3') 

            mean_k, mean_k_err = get_post_mean(post, TC, NTR, norm_const)
            logger.info('g_e_f_p4')       
        else:
            logger.info('posterior remains zero: rates set to k_bound_lo') 
            [mean_k, map_k, ci_lo, ci_hi] = [k_bound_lo for i in range(4)]
            post_grid = np.nan

        return np.asarray([mean_k, map_k, ci_lo, ci_hi]), post_grid #include post_grid for plotting post dist


    @jit(nb.float64(nb.float64[:], nb.float64[:]), nopython=True) 
    def get_chi2_jit(ntr_meas, ntr_model):
        eps = 1e-16
        chi2 = 0
        for i in range(ntr_meas.shape[0]):
            chi2 = chi2 + (ntr_model[i]-ntr_meas[i])**(2) / (ntr_model[i]+eps)
        return chi2


    @jit(nb.float64[:](types.unicode_type, types.unicode_type,
                    types.DictType(types.unicode_type, nb.float64[:]), types.unicode_type), nopython=True) 
    def get_model_pred(rr, frm, k, e_type):
        N_time = time_measured.shape[0]
        NTR_model = np.empty(N_time) 
        if e_type == '.Mean':
            idx = 0
        else:#'MAP'
            idx = 1       

        if frm == 'chr_fit':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_chr(k[rr+'T_chr'][idx], 
                                           time_measured[i])
        elif frm == 'nuc_fit':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_nuc(k[rr+'T_nuc'][idx], 
                                           time_measured[i])   
        elif frm == 'nuc_fit_from_chr':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_nuc_from_chr(k[rr+'T_chr'][idx], 
                                                    k[rr+'T_nuc'][idx], 
                                                    time_measured[i])
        elif frm == 'nuc_fit_from_nucdeg':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_nuc((k[rr+'T_nucdeg'][idx] + k[rr+'T_nucexp_from_nucdeg'][idx]), 
                                              time_measured[i])
        elif frm == 'cyto_fit_from_nucres':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_cyto(k[rr+'T_nuc'][idx], 
                                            k[rr+'T_cyto'][idx], 
                                            time_measured[i])
        elif frm == 'cyto_fit_from_chr':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_cyto_from_chr(k[rr+'T_chr'][idx], 
                                                     k[rr+'T_nucexp_from_chr'][idx], 
                                                     k[rr+'T_cyto'][idx], 
                                                     time_measured[i])
        elif frm == 'cyto_fit_from_nucdeg':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_cyto((k[rr+'T_nucdeg'][idx] + k[rr+'T_nucexp_from_nucdeg'][idx]),
                                            k[rr+'T_cyto'][idx], 
                                            time_measured[i])
        elif frm == 'poly_fit_from_nucres':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_poly(k[rr+'T_nuc'][idx], 
                                            k[rr+'T_cyto'][idx],
                                            k[rr+'T_poly_entry'][idx], 
                                            time_measured[i]) 
        elif frm == 'tot_fit':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_total_one_step(k[rr+'T_whole_cell'][idx], 
                                                      time_measured[i])               
        elif frm == 'tot_pred_from_nucres':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_total(k[rr+'T_nuc'][idx], 
                                             k[rr+'T_cyto'][idx], 
                                             time_measured[i])
        elif frm == 'tot_pred_from_chr':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_total_from_chr(k[rr+'T_chr'][idx],
                                                      k[rr+'T_nucexp_from_chr'][idx], 
                                                      k[rr+'T_cyto'][idx], 
                                                      time_measured[i])
        elif frm == 'tot_fit_from_nucdeg':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_total_from_nucdeg(k[rr+'T_nuc'][idx], 
                                                         k[rr+'T_cyto'][idx], 
                                                         k[rr+'T_nucexp_from_nucdeg'][idx],
                                                         time_measured[i])
        elif frm == 'tot_pred_from_nucdeg':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_total_from_nucdeg((k[rr+'T_nucdeg'][idx]+k[rr+'T_nucexp_from_nucdeg'][idx]), 
                                                         k[rr+'T_cyto'][idx], 
                                                         k[rr+'T_nucexp_from_nucdeg'][idx],
                                                         time_measured[i])
        return NTR_model


    logger.info('Get gene list to analyze')     
    filename = 'genes_w_rates.csv'
    genespath = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20210615')
    genes_w_rates = pd.read_csv(os.path.join(genespath, filename), 
                                             sep='\t', header=None, 
                                             squeeze=True)
    genes_w_rates = genes_w_rates[batch_size * args.start_id : min(len(genes_w_rates), batch_size * (args.start_id + 1))]


    logger.info('Run model parameter fitting for batch of genes')

    filename = 'Bayes20220228_batch'+str(args.start_id)+'.tsv'

    #initialize the fit dictionary with columns
    fits = dict()
    fits['Gene'] = []
    fits['Symbol'] = []
    for rt in RATE_TYPE:
        for rr in red_reps:
            for ts in Timescales:
                for suf in OUT_TYPES:
                    fits[rr+'.'+ts.replace('T_',rt)+suf] = []

    for r in reps:
        for fr in fracs:
            if r+fr+'top1000' in GS.keys():
                for frm in fracs_model[fr]:
                    for estimate_type in OUT_TYPES[:2]:
                        fits[red_r[r]+'.'+frm+estimate_type+'.chi2'] = []

    #Fit all genes
    count = 0
    for ensid in genes_w_rates:
        logger.info('%d %s' % (count, ensid))
        if (count % 10) == 0:
            logger.info('Write results to file %s' % filename)
            fits_df = pd.DataFrame(fits)
            fits_df.to_csv(os.path.join(outpath, filename), sep='\t',index=False)

        symbol = ''
        #initialize k_fit per gene
        k_fit = Dict.empty(key_type=types.unicode_type,
                           value_type=types.float64[:])
        for rr in red_reps:
            for ts in Timescales:
                k_fit[rr+ts] = np.asarray([np.nan for out in OUT_TYPES])


        for rr in red_reps:  
            for ts in Timescales:
                AB = dict()
                TC_TYPES_gene = dict()
                fracs_ts = Timescales2fracs[ts]
                try:
                    ABORT = False
                    for fr in fracs_ts:
                        for r in red2reps[rr]:
                            if r+fr+'top1000' in GS.keys():
                                logger.info('%s %s %s' % (rr, ts, fr))
                                TC_TYPES_gene[fr] = fit.get_tc_types_for_gene_jit(ensid, r, fr, GS.keys(), TB)
                                # The Dict.empty() constructs a typed dictionary.
                                # The key and value typed must be explicitly declared.
                                AB[fr] = Dict.empty(key_type=types.unicode_type,
                                                value_type=types.float64[:],)

                                for i in TC_TYPES_gene[fr]:
                                    tc = TC_from_jit[i]
                                    gene_idx = GS[r+fr+tc][GS[r+fr+tc]['Gene']==ensid].index[0]
                                    symbol = GS[r+fr+tc]['Symbol'][gene_idx]
                                    for ab in CI_para:
                                        gs_ab_times = [r+t+'.'+ab for t in time_id[1:]]
                                        AB[fr][tc+ab] = np.asarray(GS[r+fr+tc].loc[gene_idx, gs_ab_times], dtype='float64')
                                        logger.info('AB %s %s %s' % (tc, ab, AB[fr][tc+ab]))
                                        if np.prod(np.isfinite(AB[fr][tc+ab])) == 0:
                                            logger.warning('%s %s: abort because AB contains nonfinite data' % (r,fr))
                                            ABORT = True
                                            break
                                    if ABORT:
                                        break
                                if ABORT:
                                    break

                    if not ABORT:
                        try:
                            if ts == 'T_chr':
                                k_fit[rr + ts], _ = \
                                    get_estimates_from_post(post_chr, TC_TYPES_gene[fr], AB[fr], 1)
                            elif ts == 'T_nuc':
                                k_fit[rr + ts], _ = \
                                    get_estimates_from_post(post_nuc, TC_TYPES_gene[fr], AB[fr], 1)   
                            elif ts == 'T_nucexp_from_chr':
                                k_fit[rr + ts], _ = \
                                    get_estimates_from_post(post_nucexp_from_chr, TC_TYPES_gene, AB, 1)
                            elif ts == 'T_cyto' and \
                                sum(np.isnan(k_fit[rr + 'T_nuc'][2:])) == 0:         
                                k_fit[rr + ts], _ = \
                                    get_estimates_from_post(post_cyto, TC_TYPES_gene[fr], AB[fr], 2, 
                                                            [k_fit[rr + 'T_nuc'][2:]])
                            elif ts == 'T_poly_entry' and \
                                sum(np.isnan(k_fit[rr + 'T_nuc'][2:])) == 0 and \
                                sum(np.isnan(k_fit[rr + 'T_cyto'][2:])) == 0:
                                k_fit[rr+ ts], _ = \
                                    get_estimates_from_post(post_poly_from_nuc, TC_TYPES_gene[fr], AB[fr], 3, 
                                                                              [k_fit[rr+'T_nuc'][2:4], 
                                                                               k_fit[rr+'T_cyto'][2:4]])
                            elif ts == 'T_whole_cell':
                                k_fit[rr +ts], _ = \
                                    get_estimates_from_post(post_tot, TC_TYPES_gene[fr], AB[fr], 1)    
                            elif ts == 'T_nucexp_from_nucdeg' and \
                                sum(np.isnan(k_fit[rr + 'T_nuc'][2:])) == 0 and \
                                sum(np.isnan(k_fit[rr + 'T_cyto'][2:])) == 0:
                                k_fit[rr + ts], _ = \
                                    get_estimates_from_post(post_total_from_nucdeg, TC_TYPES_gene[fr], AB[fr], 3, 
                                                            [k_fit[rr + 'T_nuc'][2:],
                                                             k_fit[rr + 'T_cyto'][2:]])       
                            elif ts == 'T_nucdeg' and \
                                sum(np.isnan(k_fit[rr + 'T_nucexp_from_nucdeg'][2:])) == 0:
                                pnucdeg = post_nucdeg(k_fit[rr + 'T_nucexp_from_nucdeg'][2:])
                                k_fit[rr + ts], _ = \
                                    get_estimates_from_post(pnucdeg, TC_TYPES_gene[fr], AB[fr], 1) 


                        except (IndexError, KeyError, ZeroDivisionError):
                            logger.warning('except error, so nan rates')

                except (IndexError, KeyError, ZeroDivisionError):
                    logger.warning('%s %s %s: no data, so nan rates' % (rr, ts, fr))


        #append rates in fits dictionary  
        for rt in RATE_TYPE:
            for rr in red_reps:
                for ts in Timescales:
                    for i, out in enumerate(OUT_TYPES):
                        if rt == 'half_life_':
                            fits[rr+'.'+ts.replace('T_',rt)+out].append(np.log(2)*(k_fit[rr+ts][i])**(-1))
                        elif rt == 'k_':
                            fits[rr+'.'+ts.replace('T_',rt)+out].append(k_fit[rr+ts][i])
                        elif rt == 'T_':
                            fits[rr+'.'+ts+out].append((k_fit[rr+ts][i])**(-1))

        #Get predictions
        for r in reps:
            for fr in fracs:           
                if r+fr+'top1000' in GS.keys():
                    for estimate_type in OUT_TYPES[:2]:
                        gs_times = [r+t+estimate_type for t in time_id[1:]]                  
                        try:
                            TC_TYPES_gene = fit.get_tc_types_for_gene_jit(ensid, r, fr, GS.keys(), TB)

                            #get NTR timeseries as average from top1000 and bottom500
                            NTR = np.asarray([0 for _ in gs_times])
                            for i in TC_TYPES_gene:
                                tc = TC_from_jit[i]                            
                                gene_idx = GS[r+fr+tc][GS[r+fr+tc]['Gene']==ensid].index[0]
                                NTR = NTR + np.asarray(GS[r+fr+tc].loc[gene_idx, gs_times], dtype='float64')
                            NTR = NTR / len(TC_TYPES_gene)

                            for frm in fracs_model[fr]:                 
                                NTR_model = get_model_pred(red_r[r], frm, k_fit, estimate_type)
                                chi2 = get_chi2_jit(NTR, NTR_model)
                                fits[red_r[r]+'.'+frm+estimate_type+'.chi2'].append(chi2)  

                        except (IndexError, KeyError):
                            for frm in fracs_model[fr]:
                                fits[red_r[r]+'.'+frm+estimate_type+'.chi2'].append(np.nan)

        fits['Gene'].append(ensid)    
        fits['Symbol'].append(symbol)

        count+=1


    logger.info('Write results to file %s' % filename)
    fits_df = pd.DataFrame(fits)
    fits_df.to_csv(os.path.join(outpath, filename), sep='\t',index=False)
    logger.info('end')

    
if __name__ == '__main__':
    main()