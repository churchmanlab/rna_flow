#!/usr/bin/env python
# coding: utf-8

# ## Bayes Factor calculation to decide between nuc res and deg models
# Author: Robert Ietswaart  
# Date: 20220315  
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4

# Source: `Timescale_Bayes_20220222.ipynb`  
# For Subcellular Timelapse seq project: 
# Bayesian model comparison through the Bayes factor calculation to formally identify which model explains the data best: nuclear residence model or nuclear degradation.

# This script is used for batch scripts used to run for all genes together with `Bayes_Factor_20220315.sh`.   

import os
import copy
import numpy as np
import pandas as pd
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
import logging
import argparse
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
import seaborn as sns
from scipy.integrate import quad, dblquad, tplquad, cumulative_trapezoid

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
    batch_size = 25

    N_grid = 1000 #to divide [k_bound_lo,k_bound_hi] interval
    k_bound_lo = 1e-4 #1e-4 unit: min^-1: 1 per 7 days
    k_bound_hi = 1e4 #unit: min^-1: 1 per 6 ms, if too restrictive: increase to 1e6
    N_range = 8 #8 takes 3min per gene longer than 7, but 8 is same domain as calculated from posteriors.
    eps_int = 1e-3 #numerical integration precision: default: 1e-3, can go low to 1e-6

    parser = argparse.ArgumentParser(
        description='Bayes factor calculation for residence vs nucdeg model'
                    'in batches of %s genes.' % str(batch_size))
    parser.add_argument('--start_id',type=int, default=0,
                        help='row index of starting gene for batch.')
    parser.add_argument('--organism',type=str, default='m',
                        help='row index of starting gene for batch.')
    parser.add_argument('--gwr_folder',type=str, default='Bayes20210615',
                        help='row index of starting gene for batch.')   
    args = parser.parse_args()

    outpath = os.path.join('/n','groups','churchman','ri23','bseq','BayesFactor20220307')

    organisms = ['m', 'h']
    org_map = {'m': 'mouse', 'h': 'human'}
    fracs = ['nuc', 'cyto', 'tot']
    reps = ['G', 'H', 'R', 'S', 'T', 'U']
    red_reps = ['G_R', 'H_S', 'T', 'U']
    org_reps = {'m': ['G','H','R','S'], 'h': ['T', 'U']}
    org_red_reps = {'m': ['G_R','H_S'], 'h': ['T', 'U']}
    red_r = {'G': 'G_R', 'H': 'H_S', 'R': 'G_R', 'S': 'H_S', 'T': 'T', 'U': 'U'}
    red2reps = {'G_R': ['R','G'], 'H_S': ['S','H'], 'T':['T'], 'U':['U']}

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('Bayes_factor')
    log_file = os.path.join(outpath,'LogErr', 'BayesFactor20220315_%s_batch%s_py.log' % (args.organism, str(args.start_id)))
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)
    
    time_id = [str(i) for i in range(1,6)]
    time_mins = [0,15,30,60,120]
    time_measured = np.asarray(time_mins[1:])

    TC_TYPES = ['top1000','bottom500']#BM model turnover method to estimate TC
    TC_from_jit = {1: 'top1000', 0: 'bottom500'}
    CI_para = ['alpha','beta']
    OUT_TYPES = ['.Mean', '.MAP', '.0.975.quantile', '.0.025.quantile'] # 
    #the CIs must match alpha/2 and 1-alpha/2 in get_post_ci()

    Timescales = ['T_nuc',
                  'T_cyto',
                  'T_nucexp_from_nucdeg',
                  'T_nucdeg']

    GS = dict()         #GRAND-SLAM
    TB = dict()         #TopBottom genes
    B = dict()          #Bayes fits file         
    for o in organisms:  
        GS[o] = dict()   
        if o == 'm':
            path_gs = os.path.join('/n','groups','churchman','ri23','bseq','GS20210506')
            path_b = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20220222')
            filename_b = 'Bayes20220228.tsv'   

        elif o == 'h':
            path_gs = os.path.join('/n','groups','churchman','ri23','bseq','GS20210713_human')
            path_b = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20220228_human')
            filename_b = 'Bayes20220228_human.tsv'

        B[o] = pd.read_csv(os.path.join(path_b, filename_b), sep='\t')

        for r in reps:
            for fr in fracs:
                for tc in TC_TYPES:
                    filename_gs = r + '_' + fr + '_noMT_' + tc + '.csv'
                    if os.path.exists(os.path.join(path_gs, filename_gs)):
                        GS[o][r+fr+tc]= pd.read_csv(os.path.join(path_gs, filename_gs) ,index_col=0)           

                        if tc == 'top1000':
                            if o == 'm':
                                SUFFICES = ['']
                            else:
                                SUFFICES = ['genes.turnover']
                        if tc == 'bottom500':
                            SUFFICES = ['', '_below']

                        for suffix in SUFFICES:
                            filename = r + '1-' + fr + '_' + tc + suffix + '.MAPs.txt'

                            if tc == 'top1000':
                                TB[r+fr+tc]= pd.read_csv(os.path.join(path_gs, filename), 
                                                         sep='\t', header=None, 
                                                         names=['ENS_ID', 'Symbol', 'MAP'])
                            if tc == 'bottom500':
                                TB[r+fr+tc+suffix]= pd.read_csv(os.path.join(path_gs, filename), 
                                                         sep='\t', header=None, 
                                                         names=['ENS_ID', 'Symbol', 'MAP'])       



    # ### Bayes Factor calculations

    @jit(nb.float64(nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:]), nb.boolean), 
         nopython=True) 
    def pr_nuc(a, TC, NTR, pr_not_post):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_nuc(a, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            if pr_not_post: #probability for Bayes factor
                prob = prob * beta_pdf_g
            else: #posterior on rates: for increased var check
                prob = prob * beta_pdf_g * p.det_jac_nuc(a, t)
        return prob  


    @jit(nb.float64(nb.float64, nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:]), 
                    nb.boolean), nopython=True) 
    def pr_cyto(a, kD, TC, NTR, pr_not_post):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_cyto(a, kD, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            if pr_not_post: #probability for Bayes factor
                prob = prob * beta_pdf_g
            else: #posterior on rates: for increased var check
                prob = prob * beta_pdf_g * p.det_jac_cyto_from_nuc(a, kD, t)
        return prob 


    @jit(nb.float64(nb.float64, nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), 
         nopython=True) 
    def pr_tot(a, kD, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_total(a, kD, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g
            #NB: only used for Bayes factor, so no need for pr_not_post
        return prob 


    @jit(nb.float64(nb.float64, nb.float64, nb.float64, nb.int32[:], 
                    types.DictType(types.unicode_type, nb.float64[:]), nb.boolean), nopython=True) 
    def pr_total_from_nucdeg(a, kD, kE, TC, NTR, pr_not_post):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_total_from_nucdeg(a, kD, kE, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)      
            if pr_not_post: #probability for Bayes factor
                prob = prob * beta_pdf_g
            else: #posterior on rates: for increased var check
                prob = prob * beta_pdf_g * p.det_jac_total_from_nucdeg(a, kD, kE, t)
        return prob


    def get_conditional_dist(post, TC, NTR, ci, nvar, p_n_d=0):
        '''
        Calculate the marginal posterior by integration over the previous rate with 
        ci as integration boundaries.
        post: multivariate posterior
        '''
        def marg_dist(x, TC, NTR, p_n_d):
            if nvar == 2:
                return quad(post, ci[0][0], ci[0][1], args=(x, TC, NTR, p_n_d), epsabs=eps_int)[0]
            elif nvar == 3:
                return dblquad(post, ci[0][0], ci[0][1], ci[1][0], ci[1][1], args=(x, TC, NTR, p_n_d), 
                               epsabs=eps_int)[0]
        cond_dist = marg_dist
        return cond_dist


    def get_post_grid(post_, k_, TC_, NTR_):
        pg = np.empty(N_grid)
        for i in range(N_grid):
            pg[i] = post_(k_[i], TC_, NTR_, 0) #evaluate posterior in a grid to see if posterior is >0
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


    def check_NTR_using_post(post, TC_, NTR_, N_var, ci_k_prev=None):
        '''
        post : possibly multivariate posterior (function). If N_var > 1, the marginal
               distribution of the unknown rate, conditioned on ci_k_prev, will replace this
               multivariate posterior.
        N_var : number of rate arguments of the posterior.
        ci_k_prev : credible intervals of CIs of the previously fitted rates (list of tuples).
        ''' 
        if N_var > 1:
            post_new = get_conditional_dist(post, TC_, NTR_, ci_k_prev, N_var)
            post = post_new
        logger.info('g_e_f_p0')

        #Initialize the prior domain of rate
        var_scale = 1
        k_domain = np.geomspace(k_bound_lo, k_bound_hi, num=N_grid)
        post_grid, max_post = get_post_grid(post, k_domain, TC_, NTR_)
        logger.info('max post %f' % max_post)
        logger.info('g_e_f_p1')

        if max_post == 0:
            logger.info('increase variance to get nonzero posterior') 
            for vs in range(2, 5):
                NTR_scale = increase_var_NTR(TC_, NTR_, vs)
                post_scale, max_scale = get_post_grid(post, k_domain, TC_, NTR_scale)
                if max_scale > 0:
                    logger.info('max post increased: %f, var_scale %dx: %s' % (max_scale, vs, NTR_scale))
                    NTR_ = NTR_scale
                    post_grid = post_scale
                    max_post = max_scale
                    var_scale = vs
                    break
                else: 
                    logger.info('failed attempt to increase var:  var_scale %d' % vs)              
        return NTR_, var_scale


    def get_llh_nucres(TC, NTR):
        '''likelihood of nuclear residence model'''
        def integrand(k_nuc, k_cyto):
            prob = pr_nuc(k_nuc, TC['nuc'], NTR['nuc'], 1) \
                   * pr_cyto(k_nuc, k_cyto, TC['cyto'], NTR['cyto'], 1) \
                   * pr_tot(k_nuc, k_cyto, TC['tot'], NTR['tot'])
            return prob    

        #integrate over small subdomains
        llh = 0
        llh_err = 0
        for i in range(N_range):
            k_int_lo_nuc = k_bound_lo * 10**(i)
            k_int_hi_nuc = k_bound_lo * 10**(i+1)
            for j in range(N_range):
                k_int_lo_cyto = k_bound_lo * 10**(j)
                k_int_hi_cyto = k_bound_lo * 10**(j+1)
                try:
                    k_int, k_int_err = dblquad(integrand, k_int_lo_nuc, k_int_hi_nuc, 
                                                          k_int_lo_cyto, k_int_hi_cyto, epsabs=eps_int)
    #                 logger.info('llh_nucres %f %f I(_err): %f %f' % \
    #                             (k_int_lo_nuc, k_int_lo_cyto, k_int, k_int_err))
                    if np.isnan(k_int):
                        logger.warning('llh_nucres_err1 %f %f' % (k_int_lo_nuc, k_int_lo_cyto))
                        k_int = 0 
                        k_int_err = 0
                except (IndexError, KeyError, ZeroDivisionError):
                    logger.warning('llh_nucres_err2 %f %f' % (k_int_lo_nuc, k_int_lo_cyto))
                    k_int = 0 
                    k_int_err = 0
                llh += k_int
                llh_err += k_int_err
        if llh == 0:
            logger.info('llh_nucres_nb3: llh is zero') 
        llh = (k_bound_lo * 10**(N_range) - k_bound_lo)**(-2) * llh
        logger.info('llh nucres volume normalized: %s' % llh)    
        #DOES NOT WORK: domain too large
    #     llh, llh_err = dblquad(integrand, k_bound_lo, k_bound_hi, k_bound_lo, k_bound_hi, epsabs=eps_int) 
        return llh

    def get_llh_nucdeg(TC, NTR):
        '''likelihood of nuclear degradation model'''
        def integrand(k_nucdeg, k_nucexp, k_cyto):
            prob = pr_nuc((k_nucdeg + k_nucexp), TC['nuc'], NTR['nuc'], 1) \
                   * pr_cyto((k_nucdeg + k_nucexp), k_cyto, TC['cyto'], NTR['cyto'], 1) \
                   * pr_total_from_nucdeg((k_nucdeg + k_nucexp), k_cyto, k_nucexp, TC['tot'], NTR['tot'], 1)
            return prob

        #integrate over small subdomains
        llh = 0
        llh_err = 0
        for i in range(N_range):
            k_int_lo_nuc = k_bound_lo * 10**(i)
            k_int_hi_nuc = k_bound_lo * 10**(i+1)
            for j in range(N_range):
                k_int_lo_nucexp = k_bound_lo * 10**(j)
                k_int_hi_nucexp = k_bound_lo * 10**(j+1)
                for m in range(N_range):
                    k_int_lo_cyto = k_bound_lo * 10**(m)
                    k_int_hi_cyto = k_bound_lo * 10**(m+1)
                    try:
                        k_int, k_int_err = tplquad(integrand, k_int_lo_nuc, k_int_hi_nuc, 
                                                              k_int_lo_nucexp, k_int_hi_nucexp,
                                                              k_int_lo_cyto, k_int_hi_cyto,epsabs=eps_int)
    #                     logger.info('llh_nucdeg %f %f %f I(_err): %f %f' % \
    #                                 (k_int_lo_nuc, k_int_lo_cyto, k_int_lo_nucexp, k_int, k_int_err))
                        if np.isnan(k_int):
                            logger.warning('llh_nucdeg_err1 %f %f %f' % (k_int_lo_nuc, k_int_lo_nucexp, k_int_lo_cyto))
                            k_int = 0 
                            k_int_err = 0
                    except (IndexError, KeyError, ZeroDivisionError):
                        logger.warning('llh_nucdeg_err2 %f %f %f' % (k_int_lo_nuc, k_int_lo_nucexp, k_int_lo_cyto))
                        k_int = 0 
                        k_int_err = 0
                    llh += k_int
                    llh_err += k_int_err
        if llh == 0:
            logger.info('llh_nucdeg_nb3: llh is zero')       
        llh = (k_bound_lo * 10**(N_range) - k_bound_lo)**(-3) * llh
        logger.info('llh nucdeg volume normalized: %s' % llh) 

    #     #DOES NOT WORK: domain too large
    #     llh, llh_err = tplquad(integrand, k_bound_lo, k_bound_hi, 
    #                                       k_bound_lo, k_bound_hi,
    #                                       k_bound_lo, k_bound_hi, epsabs=eps_int)
    #     logger.info('llh llh_err non-volume normalized nucdeg %s %s' % (llh, llh_err))
        return llh


    def get_bayes_factor(TC, NTR, k_prev):
        '''
        TC: TC_types (top1000 and/or bottom500) for nuc, cyto and tot fractions
        NTR: alpha and betas for nuc, cyto and tot fractions
        calculate bayes_factor K:
        P(D|M0: nucres model), if nan -> K = nan, if 0 -> that is ok, output inf
        P(D|M1: nucdeg model), if nan -> K = nan
        K = P(D|M1) / P(D|M0)   
        ''' 
        out = dict()
        for fr in fracs:
            logger.info('fraction %s' % fr)
            if fr == 'nuc':
                NTR[fr], out['var_scale_' + fr] = check_NTR_using_post(pr_nuc, TC[fr], NTR[fr], 1)   
            elif fr == 'cyto' and sum(np.isnan(k_prev['nuc'])) == 0:         
                NTR[fr], out['var_scale_' + fr] = check_NTR_using_post(pr_cyto, TC[fr], NTR[fr], 2, 
                                                                       [k_prev['nuc']])
            elif fr == 'tot' and \
                sum(np.isnan(k_prev['nuc'])) == 0 and \
                sum(np.isnan(k_prev['cyto'])) == 0:
                NTR[fr], out['var_scale_' + fr] = check_NTR_using_post(pr_total_from_nucdeg, TC[fr], NTR[fr], 3, 
                                                                       [k_prev['nuc'], k_prev['cyto']])    

        out['likelihood_nucres'] = get_llh_nucres(TC, NTR)
        out['likelihood_nucdeg'] = get_llh_nucdeg(TC, NTR)
        if out['likelihood_nucres'] > 0:
            out['bayes_factor'] = out['likelihood_nucdeg'] / out['likelihood_nucres']
        else:
            if out['likelihood_nucdeg'] > 0:
                out['bayes_factor'] = np.inf
            else:# 0/0 = ill defined
                out['bayes_factor'] = np.nan 
        logger.info('bayes factor: %f' % out['bayes_factor'])               
        return out


    def get_na_out():
        out = dict()
        for fr in fracs:
            out['var_scale_' + fr] = np.nan    
        out['likelihood_nucres'] = np.nan
        out['likelihood_nucdeg'] = np.nan
        out['bayes_factor'] = np.nan
        return out
    
    

    logger.info('Get gene list to analyze')     
    filename = 'genes_w_rates.csv'
    genespath = os.path.join('/n','groups','churchman','ri23','bseq',args.gwr_folder)
    genes_w_rates = pd.read_csv(os.path.join(genespath, filename), 
                                             sep='\t', header=None, 
                                             squeeze=True)

    genes_w_rates = genes_w_rates[batch_size * args.start_id : \
                                  min(len(genes_w_rates), batch_size * (args.start_id + 1))]


    logger.info('Run bayes factor calculation for batch of genes')

    o = args.organism
    
    filename = 'Bayes_factor_20220315_'+org_map[o]+'_batch'+str(args.start_id)+'.tsv'

    #initialize the fit dictionary with columns
    fits = dict()
    fits['Gene'] = []
    fits['Symbol'] = []
    for rr in org_red_reps[o]:
        fits[rr+'.bayes_factor'] = []
        fits[rr+'.likelihood_nucres'] = []
        fits[rr+'.likelihood_nucdeg'] = []
        for fr in fracs:
            fits[rr+'.var_scale_' + fr] = []    

    #Fit all genes
    count = 0
    for ensid in genes_w_rates:
        logger.info('%d %s' % (count, ensid))
        logger.info('Write results to file %s' % filename)
        fits_df = pd.DataFrame(fits)
        fits_df.to_csv(os.path.join(outpath, filename), sep='\t',index=False)

        symbol = ''
        out = dict()
        for rr in org_red_reps[o]:  
            AB = dict()
            TC_TYPES_gene = dict()
            k_prev = dict()
            try:
                ABORT = False
                for fr in fracs:
                    for r in red2reps[rr]:
                        if r+fr+'top1000' in GS[o].keys():
                            logger.info('%s %s %s' % (r, rr, fr))
                            TC_TYPES_gene[fr] = fit.get_tc_types_for_gene_jit(ensid, r, fr, GS[o].keys(), TB)
                            logger.info('TC_types %s ' % TC_TYPES_gene[fr])
                            # The Dict.empty() constructs a typed dictionary.
                            # The key and value typed must be explicitly declared.
                            AB[fr] = Dict.empty(key_type=types.unicode_type,
                                            value_type=types.float64[:],)

                            for i in TC_TYPES_gene[fr]:
                                tc = TC_from_jit[i]
                                gene_idx = GS[o][r+fr+tc][GS[o][r+fr+tc]['Gene']==ensid].index[0]
                                symbol = GS[o][r+fr+tc]['Symbol'][gene_idx]
                                for ab in CI_para:
                                    gs_ab_times = [r+t+'.'+ab for t in time_id[1:]]
                                    AB[fr][tc+ab] = np.asarray(GS[o][r+fr+tc].loc[gene_idx, gs_ab_times], 
                                                               dtype='float64')
                                    logger.info('AB %s %s %s' % (tc, ab, AB[fr][tc+ab]))
                                    if np.prod(np.isfinite(AB[fr][tc+ab])) == 0:
                                        logger.warning('%s %s: abort because AB contains nonfinite data' % (r,fr))
                                        ABORT = True
                                        break
                                if ABORT:
                                    break
                            if ABORT:
                                break
                    if ABORT:
                        break
                    if fr in ['nuc', 'cyto']:
                        k_prev[fr] = B[o][B[o]['Gene']==ensid][[rr+'.k_'+fr+ot for ot in OUT_TYPES[2:]]].values.flatten()
                        logger.info('k_prev CI %s %s' % (fr, k_prev[fr]))
                if not ABORT:
                    try: 
                        logger.info('%s' % rr)
                        out[rr] = get_bayes_factor(TC_TYPES_gene, AB, k_prev) 
                    except (IndexError, KeyError, ZeroDivisionError):
                        logger.warning('except error, so NA bayes_factor')
                        out[rr] = get_na_out()
                else:
                    logger.warning('%s ABORT: missing data, so NA bayes_factor' % rr)
                    out[rr] = get_na_out()
            except (IndexError, KeyError, ZeroDivisionError):
                logger.warning('%s ABORT: missing data, so NA bayes_factor' % rr)
                out[rr] = get_na_out()
            logger.info('%s full output: %s' % (rr, out[rr]))

        #append output fits dictionary      
        for rr in org_red_reps[o]:
            for out_key in out[rr].keys():
                fits[rr+'.'+out_key].append(out[rr][out_key]) 
        fits['Gene'].append(ensid)    
        fits['Symbol'].append(symbol)

        count+=1

    logger.info('Write results to file %s' % filename)
    fits_df = pd.DataFrame(fits)
    fits_df.to_csv(os.path.join(outpath, filename), sep='\t',index=False)
    logger.info('end')

    
if __name__ == '__main__':
    main()