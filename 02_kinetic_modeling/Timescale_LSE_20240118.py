#!/usr/bin/env python
# coding: utf-8

# ## RNA decay fitting of models
# Author: Robert Ietswaart  
# Date: 20240118  
# License: BSD2.  
# Load modules j3dl and activate virtual environment using j4RNAdecay on O2.  
# Python v3.7.4, pytorch v1.4  
# 
# For RNA flow project: perform timescale fitting based on Grand-Slam new to total ratio MAP values for individual genes.  
# Source: `Timescale_fit_20210712_human.py`  



import os
import re
import copy
import numpy as np
import pandas as pd
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)
import math
import logging
import argparse
from scipy.optimize import least_squares
from scipy.stats import beta, pearsonr
import new_total_ratio as ntr
import fit

# from __init__ import __version__
from __init__ import default_logger_format, default_date_format

def main():
    np.random.seed(666)#12345)
    N_para = 5 
    eps = 1e-6
    k_bounds = (np.zeros(N_para)+eps,np.asarray([np.inf for i in range(N_para)]))
    k0 = abs(np.random.normal(1, 0.3, N_para))
    
    parser = argparse.ArgumentParser(
        description='Least square estimation rate fitting on subcellular timelapse seq for standard model')
    parser.add_argument('--organism',type=str, default='h',
                        help='organism: h for human or m for mouse.')
    args = parser.parse_args()
    
    o = args.organism
    
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
    
    outpath = os.path.join('/n','groups','churchman','ri23','bseq','LSE20240118')

    # Add a logger specific to the project and processing stage
    logger = logging.getLogger('LSE_rate')
    log_file = os.path.join(outpath,'LogErr', 'LSE20240118_rates_%s_py.log' % args.organism)
    formatter = logging.Formatter(default_logger_format,
                                  datefmt=default_date_format)
    log_handler = logging.FileHandler(log_file)
    log_handler.setFormatter(formatter)
    logger.addHandler(log_handler)  
    
    organisms = ['m', 'h']
    org_map = {'m': 'mouse', 'h': 'human'}
    fracs = ['chr', 'nuc', 'cyto', 'poly', 'tot']
    fracs_model = {'chr': ['chr_fit'], 
                   'nuc': ['nuc_fit'],
                   'cyto': ['cyto_fit_from_nuc'],
                   'poly': ['poly_fit_from_nuc'],
                   'tot': ['tot_fit']}
    reps = ['G', 'H', 'R', 'S', 'T', 'U']
    red_reps = ['G_R', 'H_S', 'T', 'U']
    org_reps = {'m': ['G','H','R','S'], 'h': ['T', 'U']}
    org_red_reps = {'m': ['G_R','H_S'], 'h': ['T', 'U']}
    red_r = {'G': 'G_R', 'H': 'H_S', 'R': 'G_R', 'S': 'H_S', 'T': 'T', 'U': 'U'}
    red2reps = {'G_R': ['R','G'], 'H_S': ['S','H'], 'T':['T'], 'U':['U']}
    
    time_id = [str(i) for i in range(1,6)]
    background_id = {r: '1' for r in reps}
    time_mins = [0,15,30,60,120]
    time_measured = pd.Series(time_mins[1:])
  
    TC_TYPES = ['top1000','bottom500']#BM model turnover method to estimate TC
    TC_from_jit = {1: 'top1000', 0: 'bottom500'}
    SUFFICES = {'top1000': ['genes.turnover'], 'bottom500': ['', '_below']}
    rt = 'half_life_'
    estimate_type = 'MAP'
        
    ptc_genes = dict() #List with protein-coding genes ENS IDs
    ptc_genes['h'] = os.path.join('/n', 'groups', 'churchman', 'mc348','TimelapseSeq','SeqFiles',
                                  'K562_ensGRCh38_MTmod_ptc_list.txt')
    ptc_genes['m'] = os.path.join('/n', 'groups', 'churchman', 'mc348','TimelapseSeq','SeqFiles',
                                  'NIH3T3_mm10_MTmod_ptc_list.txt')
    PTC = pd.read_csv(os.path.join(ptc_genes[o]))

    Timescales = ['T_chr',
                  'T_nuc',
                  'T_cyto',
                  'T_poly_entry',
                  'T_whole_cell']

    
    logger.info('Load GRAND-SLAM outputs')
    GS = dict()         #GRAND-SLAM
    TB = dict()         #TopBottom genes
    for r in org_reps[o]:
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


    def get_model_fit_parts(r, frm, tc):
        fit_idx = -1 
        model = None
        fixed_para = None
        if frm == 'chr_fit':
            fit_idx = 0
            model = ntr.lam_chr
        elif frm == 'nuc_fit':
            fit_idx = 1
            model = ntr.lam_nuc
        elif frm == 'cyto_fit_from_nuc':
            fit_idx = 2
            model = ntr.lam_cyto
            fixed_para = k_fit[tc+r][1]
        elif frm == 'poly_fit_from_nuc':
            fit_idx = 3
            model = ntr.lam_poly
            fixed_para = [k_fit[tc+r][1],k_fit[tc+r][2]]
        elif frm == 'tot_fit':
            fit_idx = 4
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
        elif frm == 'nuc_fit':
            model = ntr.lam_nuc
            pred_para = k_fit[tc+r][1:2]
        elif frm == 'cyto_fit_from_nuc':
            model = ntr.lam_cyto
            pred_para = [k_fit[tc+r][1], k_fit[tc+r][3],
                         k_fit[tc+r][2]]
        elif frm == 'poly_fit_from_nuc':
            model = ntr.lam_poly
            pred_para = [k_fit[tc+r][1], k_fit[tc+r][3],
                         k_fit[tc+r][2]]
        elif frm == 'tot_fit':
            model = ntr.lam_total_one_step
            pred_para = k_fit[tc+r][4:]

        return model, pred_para


    def get_chi2(ntr_meas,ntr_model):
        eps = 1e-16
        chi2 = 0
        for i in ntr_meas.index:
            chi2 = chi2 + (ntr_model[i]-ntr_meas[i])**2/(ntr_model[i]+eps)
        return chi2
 

    #Load latest genes_w_rates
    logger.info('Get gene list to analyze')  
    gpath = dict()
    gpath['h'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20231201_K562')
    gpath['m'] = os.path.join('/n','groups','churchman','ri23','bseq','Bayes20231201_3T3')
    filename = 'genes_w_rates.csv'
    genes_w_rates = pd.read_csv(os.path.join(gpath[o], filename), 
                                             sep='\t', header=None, 
                                             squeeze=True)
    
    logger.info('Initialize the fit dictionary with columns')
    fits = dict()
    fits['Gene'] = []
    fits['Symbol'] = []
    for rr in org_red_reps[o]:
        for tc in TC_TYPES:
            for ts in Timescales:
                fits[rr+'.'+tc+'.'+ts.replace('T_',rt)] = []
                
    for r in org_reps[o]:
        for tc in TC_TYPES:
            for fr in fracs:
                if r+fr+tc in GS.keys():
                    for frm in fracs_model[fr]:
                        fits[red_r[r]+'.'+tc+'.'+frm+'.chi2'] = []
                        
    logger.info('Fit LSE rates for all genes')
    filename = 'LSE_Rates_20240118_' + org_map[o]+ '.tsv'
    logger.info('Write results to file %s' % filename)
    
    count = 0
    for ensid in genes_w_rates:  
        if (count % 100) == 0:
            logger.info('%d / %d' % (count,len(genes_w_rates)))
            fits_df = pd.DataFrame(fits)
            fits_df.to_csv(os.path.join(outpath, filename), sep='\t',index=False)

        #initialize k_fit per gene
        k_fit = dict()
        NTR = dict()
        NTR_model = dict() #NTR_pred_measured
        for rr in org_red_reps[o]:
            for tc in TC_TYPES:
                k_fit[tc+rr] = np.asarray([np.nan for ts in Timescales])


        for r in org_reps[o]:
            gs_times = [r+t+' '+estimate_type for t in time_id[1:]]
            for fr in fracs:
                #check which TC conversion rate is applicable to the gene
                TC_TYPES_gene = fit.get_tc_types_for_gene_jit(ensid, r, fr, GS.keys(), TB)
                for i in TC_TYPES_gene:
                    tc = TC_from_jit[i]
                    if r+fr+tc in GS.keys():
                        try:
                            gene_idx = GS[r+fr+tc][GS[r+fr+tc]['Gene']==ensid].index[0]
                            symbol = GS[r+fr+tc]['Symbol'][gene_idx]
                            NTR[r+fr+tc] = pd.Series(list(GS[r+fr+tc].loc[gene_idx,gs_times])) 

                            for frm in fracs_model[fr]:
                                model, init_para, bounds_para, fixed_para, fit_idx = get_model_fit_parts(red_r[r], frm, tc)

                                if model is not None:
                                    try:
#                                         logger.info('ok'+r+fr+tc)
#                                         logger.info(model)
#                                         logger.info(time_measured)
#                                         logger.info(NTR[r+fr+tc])
#                                         logger.info(init_para)
#                                         logger.info(fixed_para)
#                                         logger.info(bounds_para)
#                                         logger.info(fit_idx)
                                        fit_val = least_squares(fit.calc_res, 
                                                                init_para, 
                                                                args=(model, time_measured, NTR[r+fr+tc], fixed_para),
                                                                bounds=bounds_para,
                                                                gtol=1e-14, ftol=1e-14,
                                                                loss='linear')
#                                         logger.info(fit_val.x[0])    
                                        k_fit[tc+red_r[r]][fit_idx] = fit_val.x[0]

                                        if len(TC_TYPES_gene) == 1:#apply this fitted rate to both top and bottom
                                            if tc == 'top1000':
                                                tc2 = 'bottom500'
                                            elif tc == 'bottom500':
                                                tc2 = 'top1000'
                                            k_fit[tc2+red_r[r]][fit_idx] = fit_val.x[0]
                                    except (ValueError):
#                                         logger.info('fail1'+r+fr+tc)
                                        pass
                        except (IndexError, KeyError):
#                             logger.info('fail2'+r+fr+tc)
                            pass


        #append rates in fits dictionary 
        for rr in org_red_reps[o]:
            for tc in TC_TYPES:
                for i, ts in enumerate(Timescales):
                    fits[rr+'.'+tc+'.'+ts.replace('T_',rt)].append(np.log(2)/k_fit[tc+rr][i])

        #Get predictions
        for r in org_reps[o]:
            for fr in fracs:           
                for tc in TC_TYPES:
                    if r+fr+tc in GS.keys():
                        gs_times = [r+t+' '+estimate_type for t in time_id[1:]]
                        try:
                            gene_idx = GS[r+fr+tc][GS[r+fr+tc]['Gene']==ensid].index[0]
                            symbol = GS[r+fr+tc]['Symbol'][gene_idx]
                            NTR[r+fr+tc] = pd.Series(list(GS[r+fr+tc].loc[gene_idx,gs_times])) 

                            for frm in fracs_model[fr]:                 
                                model, pred_para = get_model_pred_para(red_r[r], frm, tc)
                                NTR_model[r+frm+tc] = model(pred_para, time_measured)

                                chi2 = get_chi2(NTR[r+fr+tc], NTR_model[r+frm+tc])
                                fits[red_r[r]+'.'+tc+'.'+frm+'.chi2'].append(chi2)  
                        except (IndexError,KeyError):
                            for frm in fracs_model[fr]:
                                fits[red_r[r]+'.'+tc+'.'+frm+'.chi2'].append(np.nan)

        fits['Gene'].append(ensid)    
        fits['Symbol'].append(symbol)

        count+=1
        
    fits_df = pd.DataFrame(fits)
    fits_df.to_csv(os.path.join(outpath, filename), sep='\t',index=False)
    logger.info('end')  

if __name__ == '__main__':
    main()