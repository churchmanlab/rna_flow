#!/usr/bin/env python
# coding: utf-8

# ## Bayesian inference to obtain rate posteriors
# Author: Robert Ietswaart  
# Date: 20231201  
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4

# Source: `Timescale_Bayes_20231201.py`  
# For RNA flow project: 
# Bayesian inference to get timescale posteriors based on Grand-Slam new to total ratio posteriors for individual genes. 

# This script is used for batch scripts used to run for genes that need nucdeg rates together with `Timescale_Bayes_20240120_[K562|3T3]_iii.sh`.   


import os
import re
import copy
import numpy as np
import pandas as pd
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
import logging
import argparse
from scipy.integrate import quad, dblquad, tplquad, cumulative_trapezoid, nquad

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
    batch_size = 1 #10

    N_grid = 1000  #to divide [k_bound_lo,k_bound_hi] interval
    N_grid_var_check = 30 #to divide [k_bound_lo,k_bound_hi] interval
    grid_jiggle = np.random.normal(0, 0.01, N_grid_var_check)
    k_bound_lo = 1e-4 #1e-4 unit: min^-1: 1 per 7 days
    k_bound_hi = 1e4 #unit: min^-1: 1 per 6 ms, if too restrictive: increase to 1e6
    eps_int = 1e-3 #numerical integration precision: default: 1e-3, can go low to 1e-6
    eps_underflow = 1e-250 #minimally 1e-300. see sys
    eps_post = eps_underflow


    parser = argparse.ArgumentParser(
        description='Bayesian rate fitting on subcellular timelapse seq for standard and nucdeg model'
                    'in batches of %s genes.' % str(batch_size))
    parser.add_argument('--organism',type=str, default='h',
                        help='organism: h for human or m for mouse.')
    parser.add_argument('--oflow',type=float, default=1e-100,
                        help='overflow countering factor in integrations: 1, 1e-100, ')
    parser.add_argument('--uflow',type=float, default=1e250,
                        help='underflow countering factor in integrals: 1, 1e100, 1e250, ')
    parser.add_argument('--start_id',type=int, default=0,
                        help='row index of starting gene for batch.')
    args = parser.parse_args()
    
    o = args.organism
    oflow = args.oflow
    uflow = args.uflow
    
    path = dict() 
    path['h'] = dict()
    path['m'] = dict()
    path['h']['chr'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2021-05-16_T_U')
    path['h']['nuc'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2021-04-07_T_U')
    path['h']['cyto'] = path['h']['chr']
    path['h']['poly'] = path['h']['chr']
    path['h']['tot'] = path['h']['nuc']
    path['m']['chr'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2021-04-05_R_S')
    path['m']['nuc'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2020-10-06_G_H_fract')
    path['m']['cyto'] = path['m']['nuc']
    path['m']['poly'] = os.path.join('/n', 'groups', 'churchman', 'bms36', '2020-07-13_G_H')
    path['m']['tot'] = path['m']['poly']

    outpath = dict()
    outpath['h'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240120_K562')
    outpath['m'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20240120_3T3')

    
    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('Bayes_rate')
    log_file = os.path.join(outpath[o],'LogErr', 'Bayes20240120_rates_%s_batch%s_py.log' % 
                            (args.organism, str(args.start_id)))
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)  

    logger.info('organism %s' % o)
    logger.info('oflow %s' % oflow)
    logger.info('uflow %s' % uflow)
    
    organisms = ['m', 'h']
    org_map = {'m': 'mouse', 'h': 'human'}
    fracs = ['chr', 'nuc', 'cyto', 'poly', 'tot']
    fracs_model = {'chr': ['chr_fit'],  
                   'nuc': ['nuc_fit', 
                           'nuc_fit_from_chr_nucres', 
                           'nuc_fit_from_chr_nucdeg'],
                   'cyto': ['cyto_fit_from_nuc', 
                            'cyto_pred_from_chr_nucres', 
                            'cyto_pred_from_chr_nucdeg'],
                   'poly': ['poly_fit_from_nuc'],
                   'tot': ['tot_fit', 'tot_pred_from_nuc', 
                           'tot_pred_from_chr_nucres', 
                           'tot_pred_from_chr_nucdeg']}
    reps = ['G', 'H', 'R', 'S', 'T', 'U']
    red_reps = ['G_R', 'H_S', 'T', 'U']
    org_reps = {'m': ['G','H','R','S'], 'h': ['T', 'U']}
    org_red_reps = {'m': ['G_R','H_S'], 'h': ['T', 'U']}
    red_r = {'G': 'G_R', 'H': 'H_S', 'R': 'G_R', 'S': 'H_S', 'T': 'T', 'U': 'U'}
    red2reps = {'G_R': ['R','G'], 'H_S': ['S','H'], 'T':['T'], 'U':['U']}

    time_id = [str(i) for i in range(1,6)]
    background_id = {r: '1' for r in reps}
    time_mins = [0,15,30,60,120]
    time_measured = np.asarray(time_mins[1:])

    TC_TYPES = ['top1000','bottom500']#BM model turnover method to estimate TC
    TC_from_jit = {1: 'top1000', 0: 'bottom500'}
    CI_para = ['alpha','beta']
    OUT_TYPES = ['.Mean', '.MAP', '.0.025.quantile', '.0.975.quantile']
    SUFFICES = {'top1000': ['genes.turnover'], 'bottom500': ['', '_below']}
    #the CIs must match alpha/2 and 1-alpha/2 in get_post_ci()
    RATE_TYPES = ['half_life_','k_']
    POST_PARA_TYPES = ['var_scale', 'underflow']
    
    ptc_genes = dict() #List with protein-coding genes ENS IDs
    ptc_genes['h'] = os.path.join('/n', 'groups', 'churchman', 'mc348','TimelapseSeq','SeqFiles',
                                  'K562_ensGRCh38_MTmod_ptc_list.txt')
    ptc_genes['m'] = os.path.join('/n', 'groups', 'churchman', 'mc348','TimelapseSeq','SeqFiles',
                                  'NIH3T3_mm10_MTmod_ptc_list.txt')
    PTC = pd.read_csv(os.path.join(ptc_genes[o]))

    Timescales = ['T_chr',
                  'T_nuc',
                  'T_nucexp_from_nucres',
                  'T_nucexp_from_dist',
                  'T_cyto',
                  'T_poly_entry',
                  'T_whole_cell',
                  'T_nucdeg',
                  'T_nucexp_from_nucdeg',              
                  'T_chr_release_from_nucdeg']

    Timescales2fracs = {'T_chr': ['chr'],
                        'T_nuc': ['nuc'],
                        'T_nucexp_from_nucres': ['nuc'],
                        'T_nucexp_from_dist': ['nuc','chr'],
                        'T_cyto' : ['cyto'],
                        'T_poly_entry' : ['poly'],
                        'T_whole_cell' : ['tot'],
                        'T_nucdeg' : ['chr', 'nuc', 'cyto', 'tot'],
                        'T_nucexp_from_nucdeg': ['chr', 'nuc', 'cyto', 'tot'],
                        'T_chr_release_from_nucdeg': ['chr', 'nuc', 'cyto', 'tot']}

    logger.info('Load GRAND-SLAM outputs')
    GS = dict()         #GRAND-SLAM
    TB = dict()         #TopBottom genes
    for r in reps:
        for fr in fracs:
            for tc in TC_TYPES:
                path_gs = os.path.join(path[o][fr],'GS_2023',r+'-'+fr,re.sub(r'\d', '', tc))
                filename_gs = r + '_' + fr + '.tsv'
                if os.path.exists(os.path.join(path_gs, filename_gs)):
                    GS[r+fr+tc]= pd.read_csv(os.path.join(path_gs, filename_gs), sep='\t')
                    #filter for protein-coding genes
                    GS[r+fr+tc] = GS[r+fr+tc][GS[r+fr+tc]['Gene'].isin(PTC['Gene'])]
                    GS[r+fr+tc].sort_values(by='Gene',inplace=True, ignore_index=True)

                    for suffix in SUFFICES[tc]:
                        path_tb = os.path.join(path[o][fr],'STAR_2023',r+'1_'+fr)
                        filename_tb = r + '1_' + fr + '_' + tc + suffix + '.MAPs.txt'
                        if tc == 'top1000':
                            tb_key = r + fr + tc
                        if tc == 'bottom500':
                            tb_key = r + fr + tc + suffix
                        TB[tb_key]= pd.read_csv(os.path.join(path_tb, filename_tb), 
                                                sep='\t', header=None, 
                                                names=['ENS_ID', 'Symbol', 'MAP'])  
    

    #NB: NTR needs to be numba.typed.Dict and TC_TYPES_gene generated by calling fit.get_tc_types_for_gene_jit()
    @jit(nb.float64(nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), nopython=True) 
    def post_chr(a, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_chr(a, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_chr(a, t)
        return prob


    @jit(nb.float64(nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), nopython=True) 
    def pr_chr(a, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_chr(a, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g
        return prob  


    @jit(nb.float64(nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), nopython=True) 
    def post_nuc(kNR, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_nuc(kNR, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_nuc(kNR, t)
        return prob        


    @jit(nb.float64(nb.float64, nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), 
         nopython=True)
    def post_nucexp_from_nucres(kR, kE, TC, NTR):
        '''Relevant for non-PUNDs'''
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_nuc_from_chr(kR, kE, 0, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_nuc_from_chr_nucres(kR, kE, t)
        return prob
    
    
    @jit(nb.float64(nb.float64, nb.float64, nb.float64, nb.int32[:], 
                    types.DictType(types.unicode_type, nb.float64[:])), nopython=True)
    def pr_nuc_from_chr(a, b, kND, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_nuc_from_chr(a, b, kND, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g
        return prob

    
    def post_nucexp_from_dist(kE, TC, NTR):
        '''
        TC and NTR are now dictionaries with chr and nuc fractions as keys.
        See N8p119 for derivation.
        Relevant for non-PUNDs.
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
                    logger.warning('post_nucexp_from_dist_err1: %f %f' % (k_int, k_int_err))
                    prob = np.nan
                    prob_err = np.nan
                    break
            except (IndexError, KeyError, ZeroDivisionError):
                logger.warning('post_nucexp_from_dist_err2')
                prob = np.nan
                prob_err = np.nan
                break
            prob += k_int
            prob_err += k_int_err
        return prob

    
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


    @jit(nb.float64(nb.float64, nb.float64, nb.float64, nb.int32[:], 
                    types.DictType(types.unicode_type, nb.float64[:])), nopython=True)
    def pr_cyto_from_chr(a, b, kD, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_cyto_from_chr(a, b, kD, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g
        return prob


    @jit(nb.float64(nb.float64, nb.float64, nb.float64,
                    nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), nopython=True) 
    def post_poly(kNR, kD, kL, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_poly(kNR, kD, kL, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_poly_from_nuc(kNR, kD, kL, t)
        return prob 


    @jit(nb.float64(nb.float64, nb.int32[:], types.DictType(types.unicode_type, nb.float64[:])), nopython=True) 
    def post_tot(k_whole_cell, TC, NTR):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_total_one_step(k_whole_cell, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            prob = prob * beta_pdf_g * p.det_jac_total_one_step(k_whole_cell, t)
        return prob 


    @jit(nb.float64(nb.float64, nb.float64, nb.float64, nb.float64, nb.int32[:], 
                    types.DictType(types.unicode_type, nb.float64[:]), nb.boolean), nopython=True)
    def pr_total_from_chr(a, b, kD, kND, TC, NTR, pr_not_post):
        prob = 1
        for i in range(time_measured.shape[0]):
            t = time_measured[i]
            lam = ntr.lam_total_from_chr(a, b, kD, kND, t)
            beta_pdf_g = p.get_beta_pdf(i, lam, TC, NTR)
            if not pr_not_post: #posterior on rates
                prob = prob * beta_pdf_g * p.det_jac_total_from_chr(a, b, kD, kND, t)
            else: #probability for llh_nucres
                prob = prob * beta_pdf_g
        return prob


    @jit(nb.float64(nb.float64, nb.float64, nb.float64, nb.float64, 
                    types.DictType(types.unicode_type, nb.int32[:]), 
                    types.DictType(types.unicode_type,types.DictType(types.unicode_type, nb.float64[:])), 
                    nb.boolean, nb.boolean), nopython=True) 
    def post_joint_from_nucdeg(kR, kE, kD, kND, TC, NTR, underflow, pr_not_post=False):
        '''
        TC and NTR are now dictionaries with chr, nuc, cyto and total fractions as keys.
        See N8p202 for derivation.
        prob: joint posterior on rate vector [a, b, kD, kND] determined from above four fractions. 
        Note the order of rates as arguments is such that marginal dists and summary stats on kE and kND 
        can be determined through integration over a and kD first.
        a:= kR + kND, so da = dkR
        b:= kE + kND, so db = dkE
        '''
        if not underflow:
            prob = oflow * pr_chr((kR + kND), TC['chr'], NTR['chr']) \
                   * pr_nuc_from_chr((kR + kND), (kE + kND), kND, TC['nuc'], NTR['nuc']) \
                   * pr_cyto_from_chr((kR + kND), (kE + kND), kD, TC['cyto'], NTR['cyto']) \
                   * pr_total_from_chr((kR + kND), (kE + kND), kD, kND, TC['tot'], NTR['tot'], pr_not_post)
        else:#constant factor uflow to overcome the computer float precision of 1e-304
            prob = uflow * pr_chr((kR + kND), TC['chr'], NTR['chr']) \
                   * pr_nuc_from_chr((kR + kND), (kE + kND), kND, TC['nuc'], NTR['nuc']) \
                   * pr_cyto_from_chr((kR + kND), (kE + kND), kD, TC['cyto'], NTR['cyto']) \
                   * pr_total_from_chr((kR + kND), (kE + kND), kD, kND, TC['tot'], NTR['tot'], pr_not_post)
        return prob


    @jit(nb.float64(nb.float64, nb.float64, nb.float64, nb.float64, 
                    types.DictType(types.unicode_type, nb.int32[:]), 
                    types.DictType(types.unicode_type,types.DictType(types.unicode_type, nb.float64[:])), 
                    nb.boolean), nopython=True) 
    def post_nucdeg(kR, kE, kD, kND, TC, NTR, uf):
        '''Relevant for PUNDs'''
        return post_joint_from_nucdeg(kR, kE, kD, kND, TC, NTR, uf, False)


    @jit(nb.float64(nb.float64, nb.float64, nb.float64, nb.float64, 
                    types.DictType(types.unicode_type, nb.int32[:]), 
                    types.DictType(types.unicode_type,types.DictType(types.unicode_type, nb.float64[:])), 
                    nb.boolean), nopython=True) 
    def post_nucexp_from_nucdeg(kR, kND, kD, kE, TC, NTR, uf):
        '''Relevant for PUNDs'''
        return post_joint_from_nucdeg(kR, kE, kD, kND, TC, NTR, uf, False)

    @jit(nb.float64(nb.float64, nb.float64, nb.float64, nb.float64, 
                    types.DictType(types.unicode_type, nb.int32[:]), 
                    types.DictType(types.unicode_type,types.DictType(types.unicode_type, nb.float64[:])), 
                    nb.boolean), nopython=True) 
    def post_chr_release_from_nucdeg(kND, kE, kD, kR, TC, NTR, uf):
        '''Relevant for PUNDs'''
        return post_joint_from_nucdeg(kR, kE, kD, kND, TC, NTR, uf, False)


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
                return dblquad(post, ci[0][0], ci[0][1], ci[1][0], ci[1][1], 
                               args=(x, TC, NTR), epsabs=eps_int)[0]
        cond_dist = marg_dist
        return cond_dist


    def get_llh_nucres(TC, NTR, uf):
        '''likelihood of nuclear residence model: used to determine need for underflow correction'''   
        def optskR_nr(kE, kD, kND, TC, NTR, uf, pnp):
            return {'epsabs': eps_int, 'points': [kE, kD]}
        def optskE_nr(kD, kND, TC, NTR, uf, pnp):
            return {'epsabs': eps_int, 'points': [kD]}
        def optskD_nr(kND, TC, NTR, uf, pnp):
            return {'epsabs': eps_int} 
        #integrate over small subdomains
        N_range = 4
        llh = 0
        llh_err = 0
        for i in range(N_range):
            logger.info('llh_nucres %d / %d' % (i, N_range))
            k_R_lo_hi = [k_bound_lo * 100**(i), k_bound_lo * 100**(i+1)]
            for j in range(N_range):
                k_E_lo_hi = [k_bound_lo * 100**(j), k_bound_lo * 100**(j+1)]
                for m in range(N_range):
                    k_D_lo_hi = [k_bound_lo * 100**(m), k_bound_lo * 100**(m+1)]        
                    try:
                        k_int, k_int_err = nquad(post_joint_from_nucdeg, [k_R_lo_hi,k_E_lo_hi,k_D_lo_hi], 
                                                 args=(0, TC, NTR, uf, True), 
                                                 opts=[optskR_nr, optskE_nr, optskD_nr])#qagpe (> qagse)
                        logger.info('llh_nucres %d %d %d I(_err): %s %s' %
                                    (i, j, m, "{:e}".format(k_int), "{:e}".format(k_int_err)))
                        if np.isnan(k_int):
                            logger.warning('llh_nucres_err1 %f %f %f' % (i, j, m))
                            k_int = 0 
                            k_int_err = 0
                    except (IndexError, KeyError, ZeroDivisionError):
                        logger.warning('llh_nucres_err2 %f %f %f' % (i, j, m))
                        k_int = 0 
                        k_int_err = 0
                    llh += k_int
                    llh_err += k_int_err
                    if llh > eps_underflow:
                        logger.info('llh_nucres = %s > %s so passed this requirement with underflow %s: ' % 
                                        ("{:e}".format(llh), "{:e}".format(eps_underflow), uf))
                        return llh
        logger.info('llh nucres (not volume normalized): %s' % "{:e}".format(llh))
        return llh


    def get_post_grid(post_, TC_, NTR_, nvar, ngrid, uf=False):
        logger.info('evaluate posterior in a grid to see if posterior is > 0')
        k1 = np.geomspace(k_bound_lo, k_bound_hi, num=ngrid)
        if nvar != 4: #marginalized 1D posterior
            pg = np.empty(ngrid)
            for i in range(ngrid):
                pg[i] = post_(k1[i], TC_, NTR_)
            mp = max(pg)        
            cdf = cumulative_trapezoid(pg, k1, initial=0) 
            nc = max(cdf)
            uf = np.nan
        else: #nucdeg / nucexp_from_nucdeg / chr_release_from_nucdeg
            pg = np.empty(ngrid**nvar)
            k2 = k1 * (1 + grid_jiggle)
            k3 = k1 * (1 + 0.5 * grid_jiggle)
            k4 = k1 * (1 + 0.25 * grid_jiggle)
            for i in range(ngrid):
                for j in range(ngrid):
                    for m in range(ngrid):
                        for n in range(ngrid):                    
                            pg[((i * ngrid + j) * ngrid + m) * ngrid + n] = \
                                post_(k1[i], k2[j], k3[m], k4[n], TC_, NTR_, uf)
            mp = max(pg)
            nc = get_llh_nucres(TC_, NTR_, uf)
            if ((nc <= eps_underflow) or (mp <= eps_post)) and (not uf):
                logger.info('max post %s' % "{:e}".format(mp))
                logger.info('llh_nucres %s' % "{:e}".format(nc))           
                logger.info('possibly due to numerical underflow, '
                            're-attempt with underflow correction changed to True')
                uf = True
                for i in range(ngrid):
                    for j in range(ngrid):
                        for m in range(ngrid):
                            for n in range(ngrid):                    
                                pg[((i * ngrid + j) * ngrid + m) * ngrid + n] = \
                                    post_(k1[i], k2[j], k3[m], k4[n], TC_, NTR_, uf)
                mp = max(pg)
                nc = 1
            pg = np.nan #pg will not be used further for nucdeg/nucexp_from_nucdeg
        return mp, nc, pg, uf 


    def increase_var_NTR(TC_, NTR_, vs):
        '''Generate a new posterior by adjusting NTR alphas and betas such that 
           the variance of beta_pdf is increased by factor vs but mean of beta_pdf 
           remains unaltered: math derivation on 20220227 notes.
           NTR_: can be dict of dicts for each fraction or top/bottom single dict'''
        for fr in NTR_.keys():
            if type(NTR_[fr])==nb.typed.typeddict.Dict:
                NTR_[fr] = p.ntr_increased_var(TC_[fr], NTR_[fr], vs)
            else:
                NTR_ = p.ntr_increased_var(TC_, NTR_, vs)
                break
        return NTR_


    def check_NTR_using_post(post, TC_, NTR_, prevp, nvar, ngrid):
        '''
        post : posterior (function).
        prevp: previously estimated post parameters (dict with possible keys: var_scale, underflow)
        nvar: number of rate arguments of the posterior. nvar == 4 <=> nucdeg or nucexp_from_nucdeg with 
        multivariate posterior. In all other cases post is 1 dimensional.
        ''' 
        if nvar != 4:#all non-nucdeg model rates
            if (prevp is not None) and ('var_scale' in prevp.keys()): 
                var_scale = prevp['var_scale']
                logger.info('var_scale %d' % var_scale)
                if var_scale > 1:
                    NTR_ = increase_var_NTR(TC_, NTR_, var_scale)
            else:
                prevp = dict()
                var_scale = 1 

            max_post, norm_const, post_grid, uf = get_post_grid(post, TC_, NTR_, nvar, ngrid)
            logger.info('g_e_f_p1')
            logger.info('max post %s' % "{:e}".format(max_post))
            logger.info('norm_constant %s' % "{:e}".format(norm_const))

            if (max_post <= eps_post) or (norm_const <= eps_underflow):
                logger.info('increase variance to get nonzero posterior') 
                for vs in range((var_scale + 1), 5):
                    NTR_scale = increase_var_NTR(TC_, NTR_, vs)
                    max_scale, norm_scale, post_grid, uf = get_post_grid(post, TC_, NTR_scale, nvar, ngrid)
                    if (max_scale > eps_post) and (norm_scale > eps_underflow):
                        logger.info('max post increased: %s, var_scale %dx norm_const %s: %s' % 
                                    ("{:e}".format(max_scale), vs, "{:e}".format(norm_scale), NTR_scale))
                        NTR_ = NTR_scale
                        max_post = max_scale
                        norm_const = norm_scale
                        var_scale = vs
                        break
                    else: 
                        logger.info('failed attempt to increase var: var_scale %d' % vs)
            prevp['var_scale'] = var_scale

        else:#nucdeg / nucexp_from_nucdeg / chr_release_from_nucdeg
            for ts in ['T_chr', 'T_nuc', 'T_cyto', 'T_whole_cell']:
                fr = Timescales2fracs[ts][0]
                var_scale = prevp[ts]['var_scale']
                logger.info('ts %s var_scale %d' % (ts, var_scale) )
                if var_scale > 1:
                    NTR_[fr] = increase_var_NTR(TC_[fr], NTR_[fr], var_scale)
            if 'T_nucdeg' in prevp.keys():
                logger.info('g_e_f_p1: NTR and underflow already determined in nucdeg previously')
                uf = prevp['T_nucdeg']['underflow']
                logger.info('underflow %s' % uf)
                max_post = 1
                norm_const = 1 
                post_grid = np.nan
            else:
                max_post, norm_const, post_grid, uf = get_post_grid(post, TC_, NTR_, nvar, ngrid) 
                logger.info('g_e_f_p1 nucdeg')
                logger.info('max post %s' % "{:e}".format(max_post))
                logger.info('llh_nucres %s' % "{:e}".format(norm_const))
                logger.info('underflow %s' % uf)

                if (max_post <= eps_post) or (norm_const <= eps_underflow):
                    logger.warning('Multivar posterior <= eps_post: indicates data is internally inconsistent.' 
                                   'This should not be happening, investigate gene manually.')
            prevp['var_scale'] = np.nan
        prevp['underflow'] = uf

        return NTR_, prevp, max_post, norm_const, post_grid


    def get_post_ci(k_domain, post_grid, norm_const):
        alpha = 0.05
        #CDF via cumulative trapezoid
        cdf = cumulative_trapezoid(post_grid, k_domain, initial=0) 
        cdf = cdf / norm_const
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
        return ci_lo, ci_hi


    def get_post_from_nucdeg(post, TC, NTR, uf, ci):
        if post == post_nucdeg:
            def optskR(kE, kD, kND, TC, NTR, uf):
                return {'epsabs': eps_int, 'points': [kE, (kD - kND)]}
            def optskE(kD, kND, TC, NTR, uf):
                return {'epsabs': eps_int, 'points': [(kD - kND)]}
            def optskD(kND, TC, NTR, uf):
                return {'epsabs': eps_int}
            opts_integrand = [optskR, optskE, optskD]
            def kR_lo_hi(kE, kD, kND, TC, NTR, uf):
                '''ci[0][0] = a_lo
                   ci[0][1] = a_hi'''
                return [max(ci[0][0] - kND, k_bound_lo), max(ci[0][1] - kND, k_bound_lo)]  
        elif post == post_nucexp_from_nucdeg:
            def optskR(kND, kD, kE, TC, NTR, uf):
                return {'epsabs': eps_int, 'points': [kE, (kD - kND)]}
            def optskND(kD, kE, TC, NTR, uf):
                return {'epsabs': eps_int, 'points': [(kD - kE)]}
            def optskD(kE, TC, NTR, uf):
                return {'epsabs': eps_int}
            opts_integrand = [optskR, optskND, optskD]
            def kR_lo_hi(kND, kD, kE, TC, NTR, uf):
                '''ci[0][0] = a_lo
                   ci[0][1] = a_hi'''
                return [max(ci[0][0] - kND, k_bound_lo), max(ci[0][1] - kND, k_bound_lo)] 
        elif post == post_chr_release_from_nucdeg:
            def optskND(kE, kD, kR, TC, NTR, uf):
                return {'epsabs': eps_int, 'points': [(kD - kE),(kD - kR)]}
            def optskE(kD, kR, TC, NTR, uf):
                return {'epsabs': eps_int, 'points': [kR]}
            def optskD(kR, TC, NTR, uf):
                return {'epsabs': eps_int}  
            opts_integrand = [optskND, optskE, optskD]  
            def kND_lo_hi(kE, kD, kR, TC, NTR, uf):
                '''ci[0][0] = a_lo
                   ci[0][1] = a_hi'''
                return [max(ci[0][0] - kR, k_bound_lo), max(ci[0][1] - kR, k_bound_lo)] 
        else:
            logger.warning('post not identified, return NA')
            return np.nan, np.nan, np.nan, np.nan

        kD_lo_hi = ci[1]

        #initialize domain and posterior
        N_grid_nucdeg = 100
        k_domain = np.geomspace(k_bound_lo, k_bound_hi, num=N_grid_nucdeg)
        post_grid = np.empty(N_grid_nucdeg)

        #integrate over small subdomains and sum   
        N_range = 4 #does not make a difference with 8 and has lower runtime
        for i in range(N_grid_nucdeg):
            logger.info('post_from_nucdeg %d / %d' % (i, N_grid_nucdeg))
            prob = 0  
            prob_err = 0     
            for j in range(N_range):
                k2_lo_hi = [k_bound_lo * 100**(j), k_bound_lo * 100**(j+1)] #kR, kE or kND depend on post
                try:
                    if post == post_nucdeg:
                        bounds_integrand = [kR_lo_hi, k2_lo_hi, kD_lo_hi]
                    elif post == post_nucexp_from_nucdeg: 
                        bounds_integrand = [kR_lo_hi, k2_lo_hi, kD_lo_hi]
                    elif post == post_chr_release_from_nucdeg:
                        bounds_integrand = [kND_lo_hi, k2_lo_hi, kD_lo_hi]                
                    p_int, p_int_err = nquad(post, bounds_integrand,
                                             args=(k_domain[i], TC, NTR, uf), 
                                             opts=opts_integrand)#qagpe
                    logger.info('p_int %s %s' % ("{:e}".format(p_int),"{:e}".format(p_int_err)))
                    if np.isnan(p_int):
                        logger.warning('prob_err1 %f %f' % (i, j))
                        prob = np.nan
                        prob_err = np.nan
                        break
                    prob += p_int
                    prob_err += p_int_err
                except (IndexError, KeyError, ZeroDivisionError):
                    logger.warning('prob_err2 %f %f' % (i, j))
                    prob = np.nan
                    prob_err = np.nan
                    break
            if prob < 0:
                logger.warning('prob3: prob is smaller than 0: set to 0')
                prob = 0
                prob_err = np.nan
            post_grid[i] = prob

        cdf = cumulative_trapezoid(post_grid, k_domain, initial=0) 
        nc = max(cdf)
        mp = max(post_grid)
        logger.info('norm_constant from_nucdeg %s' % "{:e}".format(nc))
        logger.info('max_post from_nucdeg %s' % "{:e}".format(mp))
        return post_grid, nc, mp


    def get_post_mean(post, TC, NTR, norm_const):
        def integrand(x):
            return x * post(x, TC, NTR) / norm_const

        #integrate over small subdomains
        N_range = 8
        k_mean = 0
        k_mean_err = 0
        for i in range(N_range):
            k_int_lo = k_bound_lo * 10**(i)
            k_int_hi = k_bound_lo * 10**(i+1)
            try:
                k_int, k_int_err = quad(integrand, k_int_lo, k_int_hi, epsabs=eps_int)
                logger.info('mean %s %s' % ("{:e}".format(k_int),"{:e}".format(k_int_err)))
                if np.isnan(k_int):
                    logger.warning('mean_err1 %f %f' % (k_int_lo, k_int_hi))
                    k_mean = np.nan
                    k_mean_err = np.nan
                    break
                k_mean += k_int
                k_mean_err += k_int_err
            except (IndexError, KeyError, ZeroDivisionError):
                logger.warning('mean_err2 %f %f' % (k_int_lo, k_int_hi))
                k_mean = np.nan
                k_mean_err = np.nan
                break
        if k_mean < k_bound_lo:
            logger.warning('mean_err3: mean < k_bound_lo, set to NA')
            k_mean = np.nan
            k_mean_err = np.nan       
        if k_mean > k_bound_hi:
            logger.warning('mean_err4: mean > k_bound_hi, set to NA')
            k_mean = np.nan
            k_mean_err = np.nan     
        return k_mean, k_mean_err


    def get_estimates_from_post(post, TC, NTR, N_var, ci_k_prev=None, prevp=None):
        '''
        post : possibly multivariate posterior (function). If N_var > 1, the marginal
               distribution of the unknown rate, conditioned on ci_k_prev, will replace this
               multivariate posterior. 
        N_var : number of rate arguments of the posterior. 
        N_var == 4 implies nucdeg or nucexp_from_nucdeg: rates are not separable from 
        compartmental posteriors, unlike all other cases.
        prevp: previously estimated post parameters (dict[keys: var_scale, underflow] or 
        dict[keys:timescales] of dicts[keys: var_scale, underflow])
        ci_k_prev : credible intervals of CIs of the previously fitted rates (list of tuples).
        '''
        [mean_k, map_k, ci_lo, ci_hi] = [np.nan for i in range(4)]
        if N_var != 4:
            if N_var > 1:
                post_new = get_conditional_dist(post, TC, NTR, ci_k_prev, N_var)
                post = post_new
            logger.info('g_e_f_p0')

            NTR, postp, max_post, norm_const, post_grid = check_NTR_using_post(post, TC, NTR, prevp, 
                                                                               N_var, N_grid)

            if (max_post > eps_post) and (norm_const > eps_underflow) and np.isfinite(norm_const):

                k_domain = np.geomspace(k_bound_lo, k_bound_hi, num=N_grid)
                map_k = k_domain[np.argmax(post_grid)]
                logger.info('g_e_f_p2') 

                ci_lo, ci_hi = get_post_ci(k_domain, post_grid, norm_const)
                logger.info('g_e_f_p3') 

                mean_k, mean_k_err = get_post_mean(post, TC, NTR, norm_const)
                if np.isnan(mean_k):
                    logger.info('Mean = NA: could not be numerically calculated reliably.'
                                'Set mean equal to MAP,'
                                'since MAP is remaining best estimate for mean.')
                    mean_k = map_k
                logger.info('g_e_f_p4')       

        else:
            logger.info('g_e_f_p0: nucdeg / nucexp_from_nucdeg / chr_release_from_nucdeg')

            NTR, postp, max_post, llh_nucres, _ = check_NTR_using_post(post, TC, NTR, prevp, 
                                                                       N_var, N_grid_var_check) 

            if (max_post > eps_post) and (llh_nucres > eps_underflow):

                post_grid, norm_const, max_post = get_post_from_nucdeg(post, TC, NTR, 
                                                                       postp['underflow'], ci_k_prev)            

                if (max_post > eps_post) and (norm_const > eps_underflow) and np.isfinite(norm_const): 
                    k_domain = np.geomspace(k_bound_lo, k_bound_hi, num=post_grid.shape[0])
                    map_k = k_domain[np.argmax(post_grid)]
                    logger.info('g_e_f_p2') 

                    ci_lo, ci_hi = get_post_ci(k_domain, post_grid, norm_const)
                    logger.info('g_e_f_p3') 

                    logger.info('Mean = NA: could not be calculated reliably. '
                                'Set mean equal to MAP,'
                                'since MAP is remaining best estimate for mean.')
                    mean_k = map_k
                    logger.info('g_e_f_p4')

        if max_post <= eps_post:
            logger.info('posterior remains zero: rates set to NA') 
            [mean_k, map_k, ci_lo, ci_hi] = [np.nan for i in range(4)]

        return np.asarray([mean_k, map_k, ci_lo, ci_hi]), postp

    
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
                NTR_model[i] = ntr.lam_chr(k[rr+'T_chr'][idx], time_measured[i])
        elif frm == 'nuc_fit':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_nuc(k[rr+'T_nuc'][idx], time_measured[i])   
        elif frm == 'nuc_fit_from_chr_nucres':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_nuc_from_chr(k[rr+'T_chr'][idx], k[rr+'T_nucexp_from_nucres'][idx], 
                                                    0, time_measured[i])
        elif frm == 'nuc_fit_from_chr_nucdeg':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_nuc_from_chr((k[rr+'T_chr_release_from_nucdeg'][idx]+k[rr+'T_nucdeg'][idx]), 
                                                    (k[rr+'T_nucexp_from_nucdeg'][idx]+k[rr+'T_nucdeg'][idx]),
                                                    k[rr+'T_nucdeg'][idx], time_measured[i])
        elif frm == 'cyto_fit_from_nuc':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_cyto(k[rr+'T_nuc'][idx], k[rr+'T_cyto'][idx], time_measured[i])
        elif frm == 'cyto_pred_from_chr_nucres':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_cyto_from_chr(k[rr+'T_chr'][idx], k[rr+'T_nucexp_from_nucres'][idx], 
                                                     k[rr+'T_cyto'][idx], time_measured[i])
        elif frm == 'cyto_pred_from_chr_nucdeg':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_cyto_from_chr((k[rr+'T_chr_release_from_nucdeg'][idx]+k[rr+'T_nucdeg'][idx]), 
                                                     (k[rr+'T_nucexp_from_nucdeg'][idx]+k[rr+'T_nucdeg'][idx]),
                                                     k[rr+'T_cyto'][idx], time_measured[i])
        elif frm == 'poly_fit_from_nuc':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_poly(k[rr+'T_nuc'][idx], k[rr+'T_cyto'][idx],
                                            k[rr+'T_poly_entry'][idx], time_measured[i]) 
        elif frm == 'tot_fit':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_total_one_step(k[rr+'T_whole_cell'][idx], time_measured[i])  
        elif frm == 'tot_pred_from_nuc':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_total(k[rr+'T_nuc'][idx], k[rr+'T_cyto'][idx], time_measured[i])
        elif frm == 'tot_pred_from_chr_nucres':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_total_from_chr(k[rr+'T_chr'][idx], k[rr+'T_nucexp_from_nucres'][idx], 
                                                      k[rr+'T_cyto'][idx], 0, time_measured[i])
        elif frm == 'tot_pred_from_chr_nucdeg':
            for i in range(N_time):
                NTR_model[i] = ntr.lam_total_from_chr((k[rr+'T_chr_release_from_nucdeg'][idx]+k[rr+'T_nucdeg'][idx]),
                                                      (k[rr+'T_nucexp_from_nucdeg'][idx]+k[rr+'T_nucdeg'][idx]), 
                                                      k[rr+'T_cyto'][idx], k[rr+'T_nucdeg'][idx], time_measured[i])
        return NTR_model  
    
  

    logger.info('Get gene list to analyze')     
    filename = 'genes_w_nucdeg_reruns.tsv' #'genes_w_rates.csv'
    genes_w_rates = pd.read_csv(os.path.join(outpath[o], filename), 
                                             sep='\t', header=None, 
                                             squeeze=True)

    genes_w_rates = genes_w_rates[batch_size * args.start_id : \
                                  min(len(genes_w_rates), batch_size * (args.start_id + 1))] #run i
#     genes_w_rates = genes_w_rates[args.start_id : min(len(genes_w_rates), args.start_id + args.n_genes)] #run >= ii


    logger.info('Run model parameter fitting for batch of genes')
    
    filename = 'Bayes_Rates_20240120_'+org_map[o]+'_batch'+str(args.start_id)+'.tsv' #run i
#     filename = 'Bayes_Rates_20240120_'+org_map[o]+'_batch'+str(args.suffix_id)+'.tsv' #run >= ii
    
    #initialize the fit dictionary with columns
    fits = dict()
    fits['Gene'] = []
    fits['Symbol'] = []
    for rt in RATE_TYPES:
        for rr in org_red_reps[o]:
            for ts in Timescales:
                for suf in OUT_TYPES:
                    fits[rr+'.'+ts.replace('T_',rt)+suf] = []
                for suf in POST_PARA_TYPES:
                    fits[rr+'.'+ts.replace('T_',rt)+'.'+suf] = []

    for r in org_reps[o]:
        for fr in fracs:
            if r+fr+'top1000' in GS.keys():
                for frm in fracs_model[fr]:
                    for estimate_type in OUT_TYPES[:2]:
                        fits[red_r[r]+'.'+frm+estimate_type+'.chi2'] = []


    #Fit all genes
    count = 0
    for ensid in genes_w_rates:        
        logger.info('%d %s' % (count, ensid))
        logger.info('Write results to file %s' % filename)
        fits_df = pd.DataFrame(fits)
        fits_df.to_csv(os.path.join(outpath[o], filename), sep='\t',index=False)

        symbol = ''
        #initialize k_fit per gene
        k_fit = Dict.empty(key_type=types.unicode_type,
                           value_type=types.float64[:])
        ppara = dict()

        for rr in org_red_reps[o]:
            for ts in Timescales:
                k_fit[rr+ts] = np.asarray([np.nan for out in OUT_TYPES])
                ppara[rr+ts] = {out: np.nan for out in POST_PARA_TYPES}

        for rr in org_red_reps[o]:  
            for ts in Timescales:
                AB = Dict.empty(key_type=types.unicode_type,
                                value_type=types.DictType(types.unicode_type, types.float64[:]))
                TC_TYPES_gene = Dict.empty(key_type=types.unicode_type,
                                           value_type=types.int32[:])
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
                                        gs_ab_times = [r+t+' '+ab for t in time_id[1:]]
                                        AB[fr][tc+ab] = np.asarray(GS[r+fr+tc].loc[gene_idx, gs_ab_times], 
                                                                   dtype='float64')
                                        logger.info('AB %s %s %s' % (tc, ab, AB[fr][tc+ab]))
                                        if np.prod(np.isfinite(AB[fr][tc+ab])) == 0:
                                            logger.warning('%s %s: abort because AB contains nonfinite data' % 
                                                           (r,fr))
                                            ABORT = True
                                            break
                                    if ABORT:
                                        break
                                if ABORT:
                                    break

                    if not ABORT:
                        try:
                            if ts == 'T_chr':
                                k_fit[rr + ts], ppara[rr + ts] = \
                                    get_estimates_from_post(post_chr, TC_TYPES_gene['chr'], AB['chr'], 1)
                            elif ts == 'T_nuc':
                                k_fit[rr + ts], ppara[rr + ts] = \
                                    get_estimates_from_post(post_nuc, TC_TYPES_gene['nuc'], AB['nuc'], 1)   
#                             elif ts == 'T_nucexp_from_nucres' and \
#                                 sum(np.isnan(k_fit[rr + 'T_chr'][2:])) == 0:
#                                 k_fit[rr + ts], ppara[rr + ts] = \
#                                     get_estimates_from_post(post_nucexp_from_nucres, TC_TYPES_gene['nuc'], 
#                                                             AB['nuc'], 2, [k_fit[rr + 'T_chr'][2:]], 
#                                                             ppara[rr + 'T_nuc']) 
#                             elif ts == 'T_nucexp_from_dist':
#                                 k_fit[rr + ts], _ = \
#                                     get_estimates_from_post(post_nucexp_from_dist, TC_TYPES_gene, AB, 1)
                            elif ts == 'T_cyto' and \
                                sum(np.isnan(k_fit[rr + 'T_nuc'][2:])) == 0:         
                                k_fit[rr + ts], ppara[rr + ts] = \
                                    get_estimates_from_post(post_cyto, TC_TYPES_gene['cyto'], AB['cyto'], 2, 
                                                            [k_fit[rr + 'T_nuc'][2:]])
#                             elif ts == 'T_poly_entry' and \
#                                 sum(np.isnan(k_fit[rr + 'T_nuc'][2:])) == 0 and \
#                                 sum(np.isnan(k_fit[rr + 'T_cyto'][2:])) == 0:
#                                 k_fit[rr + ts], ppara[rr + ts] = \
#                                     get_estimates_from_post(post_poly, TC_TYPES_gene['poly'], AB['poly'], 3, 
#                                                                               [k_fit[rr+'T_nuc'][2:], 
#                                                                                k_fit[rr+'T_cyto'][2:]])    
                            elif ts == 'T_whole_cell':
                                k_fit[rr + ts], ppara[rr + ts] = \
                                    get_estimates_from_post(post_tot, TC_TYPES_gene['tot'], AB['tot'], 1)  
                            elif ts == 'T_nucdeg'  and \
                                sum(np.isnan(k_fit[rr + 'T_chr'][2:])) == 0 and \
                                sum(np.isnan(k_fit[rr + 'T_cyto'][2:])) == 0:
                                k_fit[rr + ts], ppara[rr + ts] = \
                                    get_estimates_from_post(post_nucdeg, TC_TYPES_gene, AB, 4,
                                                            [k_fit[rr+'T_chr'][2:],
                                                             k_fit[rr+'T_cyto'][2:]],
                                                            {ts: ppara[rr + ts] \
                                                             for ts in ['T_chr','T_nuc', 
                                                                        'T_cyto', 'T_whole_cell']}) 
                            elif ts == 'T_nucexp_from_nucdeg' and \
                                sum(np.isnan(k_fit[rr + 'T_chr'][2:])) == 0 and \
                                sum(np.isnan(k_fit[rr + 'T_cyto'][2:])) == 0:
                                k_fit[rr + ts], ppara[rr + ts] = \
                                    get_estimates_from_post(post_nucexp_from_nucdeg, TC_TYPES_gene, AB, 4, 
                                                            [k_fit[rr + 'T_chr'][2:],
                                                             k_fit[rr + 'T_cyto'][2:]],
                                    {ts: ppara[rr + ts] for ts in ['T_chr', 'T_nuc', 'T_cyto', 
                                                                   'T_whole_cell', 'T_nucdeg']})     
                            elif ts == 'T_chr_release_from_nucdeg' and \
                                sum(np.isnan(k_fit[rr + 'T_chr'][2:])) == 0 and \
                                sum(np.isnan(k_fit[rr + 'T_cyto'][2:])) == 0:
                                k_fit[rr + ts], ppara[rr + ts] = \
                                    get_estimates_from_post(post_chr_release_from_nucdeg, TC_TYPES_gene, 
                                                            AB, 4, [k_fit[rr + 'T_chr'][2:],
                                                                    k_fit[rr + 'T_cyto'][2:]],
                                    {ts: ppara[rr + ts] for ts in ['T_chr', 'T_nuc', 'T_cyto', 
                                                                   'T_whole_cell', 'T_nucdeg']})   
                            else:
                                logger.warning('%s does not execute' % ts)
                        except (IndexError, KeyError, ZeroDivisionError):
                            logger.warning('%s %s %s: except error: ABORT = True, so nan rates' % (rr, ts, fr))
                except (IndexError, KeyError, ZeroDivisionError):
                    logger.warning('%s %s: no AB data, so nan rates' % (rr, ts))


        #append rates in fits dictionary  
        for rt in RATE_TYPES:
            for rr in org_red_reps[o]:
                for ts in Timescales:
                    for i, out in enumerate(OUT_TYPES):
                        if rt == 'half_life_':
                            if i < 2:#Mean and MAP
                                fits[rr+'.'+ts.replace('T_',rt)+out].append(np.log(2)*(k_fit[rr+ts][i])**(-1))
                            elif i == 2:#CI lo of half-life: use  CI hi of rate, because of inverse calculation
                                fits[rr+'.'+ts.replace('T_',rt)+out].append(np.log(2)*(k_fit[rr+ts][3])**(-1))
                            else:#CI hi of half-life: use  CI lo of rate, because of inverse calculation
                                fits[rr+'.'+ts.replace('T_',rt)+out].append(np.log(2)*(k_fit[rr+ts][2])**(-1))
                        elif rt == 'k_':
                            fits[rr+'.'+ts.replace('T_',rt)+out].append(k_fit[rr+ts][i])
                    for out in POST_PARA_TYPES:
                        fits[rr+'.'+ts.replace('T_',rt)+'.'+out].append(ppara[rr+ts][out])

#         #Get predictions
#         for r in org_reps[o]:
#             for fr in fracs:           
#                 if r+fr+'top1000' in GS.keys():
#                     for estimate_type in OUT_TYPES[:2]:
#                         gs_times = [r+t+estimate_type.replace('.',' ') for t in time_id[1:]]                  
#                         try:
#                             TC_TYPES_gene = fit.get_tc_types_for_gene_jit(ensid, r, fr, GS.keys(), TB)

#                             #get NTR timeseries as average from top1000 and bottom500
#                             NTR = np.asarray([0 for _ in gs_times])
#                             for i in TC_TYPES_gene:
#                                 tc = TC_from_jit[i]                            
#                                 gene_idx = GS[r+fr+tc][GS[r+fr+tc]['Gene']==ensid].index[0]
#                                 NTR = NTR + np.asarray(GS[r+fr+tc].loc[gene_idx, gs_times], dtype='float64')
#                             NTR = NTR / len(TC_TYPES_gene)

#                             for frm in fracs_model[fr]:                 
#                                 NTR_model = get_model_pred(red_r[r], frm, k_fit, estimate_type)
#                                 chi2 = get_chi2_jit(NTR, NTR_model)
#                                 fits[red_r[r]+'.'+frm+estimate_type+'.chi2'].append(chi2)  
#                         except (IndexError, KeyError):
#                             for frm in fracs_model[fr]:
#                                 fits[red_r[r]+'.'+frm+estimate_type+'.chi2'].append(np.nan)

        fits['Gene'].append(ensid)    
        fits['Symbol'].append(symbol)
        count+=1

    logger.info('Write results to file %s' % filename)
    fits_df = pd.DataFrame(fits)
    fits_df.to_csv(os.path.join(outpath[o], filename), sep='\t',index=False)
    logger.info('end')   
    
    
if __name__ == '__main__':
    main()
