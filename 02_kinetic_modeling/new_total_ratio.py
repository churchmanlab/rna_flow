#!/usr/bin/env python
# coding: utf-8

# ## new_total_ratio.py
# Author: Robert Ietswaart
# Date: 20201201
# License: BSD2.  
# Python v3.7.4 
# 
# For Brendan's project: fit time scales to his subcellular timelapse seq using observed new to total ratio (Lambda)
# of the various subcellular RNA fractions
# See Notebook N7p51-56,p174-179,p181-191


import numpy as np
from scipy.stats import binom
import lmfit

def lam_chr(k, t, k_para=None):
    """ Fit based on chr compartment
        t: time
        k: rates list [kR, kE, kL, kD]
        kR: RNA Release rate from chromatin into nucleoplasm
        kch: RNA degradation rate at chromatin, set to zero
    """
    kR = k[0]
    kch = 0
    Lambda = 1 - np.exp(-(kR + kch) * t)
    return Lambda

# def lam_nucpl(k, t, k_para=None):
#     """ Fit based on chr and nucpl compartments
#         t: time
#         k: rates list [kR, kE, kL, kD]
#         kR: RNA Release rate from chromatin into nucleoplasm
#         kE: RNA export rate from nucleoplasm into cytoplasm
#         kch: RNA degradation rate at chromatin, set to zero
#         knp: RNA degradation rate in nucleoplasm, set to zero
#     """
#     if k_para == None:
#         kR = k[0]
#         kE = k[1]
#     else:
#         kR = k_para
#         kE = k[0]
#     kch = 0
#     knp = 0
#     a = kch + kR
#     b = knp + kE
#     if a == b:#prevent divide by zero
#         Lambda = 1 - np.exp(-a * t)
#     else:
#         Phi = a / (a - b)
#         Omega = -b / (a - b) #1-Phi
#         Lambda = Phi * (1 - np.exp(-b * t)) + \
#                  Omega * (1 - np.exp(-a * t))
#     return Lambda

def lam_nuc(k, t, k_para=None):
    """ Fit based on only nuc compartment data.
        t: time
        k: rates list [kE, kL, kD]
        kE: RNA export rate from nucleus into cytoplasm
    """
    kE = k[0]
    Lambda = 1 - np.exp(-kE * t)
    return Lambda

# def lam_nuc_from_chr(k, t, k_para=None):
#     """ Prediction as sum of chr and nucpl compartments
#         t: time
#         k: rates list [kR, kE, kL, kD]
#         kR: RNA Release rate from chromatin into nucleoplasm
#         kE: RNA export rate from nucleoplasm into cytoplasm
#         kch: RNA degradation rate at chromatin, set to zero
#         knp: RNA degradation rate in nucleoplasm, set to zero
#     """
#     if k_para == None:
#         kR = k[0]
#         kE = k[1]
#     else:
#         kR = k_para
#         kE = k[0]
#     kch = 0
#     knp = 0
#     a = kch + kR
#     b = knp + kE
#     if a == b:
#         Lambda = 1 - np.exp(-a * t)
#     else:
#         prefactor = ((kR + b) * (a - b))**(-1)
#         Phi = -b * (kR + b - a)
#         Omega = kR * a
#         Lambda = prefactor * ( Phi * (1 - np.exp(-a * t)) + \
#                  Omega * (1 - np.exp(-b * t)) )
#     return Lambda

def lam_cyto(k, t, k_para=None):
    """ Fit based on nuc, cyto compartment data.
        t: time
        k: rates list [kE, kL, kD] or [kD] when fitting kD only
        in that case k_para == kE
        kE: RNA Export rate from nucleus into cytoplasm
        kD: (cytoplasmic messenger) RNA degradation rate
    """
    if k_para == None:
        if len(k) == 2:#for Bayes
            kE = k[0]
            kD = k[1]            
        else:#maintain backwards compatibility
            kE = k[0]
            kD = k[2]
    else:
        kE = k_para
        kD = k[0]
    if kD == kE:
        Lambda = 1 - np.exp(-kD * t)
    else:
        Phi = kD / (kD - kE)
        Omega = -kE / (kD - kE) #1-Phi
        Lambda = Phi * (1 - np.exp(-kE * t)) + \
                 Omega * (1 - np.exp(-kD * t))
    return Lambda

# def lam_cyto_from_chr(k, t, k_para=None):
#     """ Fit based on chr, nucpl and cyto compartments
#         t: time
#         k: rates list [kR, kE, kL, kD] or [kD] when fitting kD only
#         in that case k_para == [kR, kE]
#         kR: RNA Release rate from chromatin into nucleoplasm
#         kE: RNA export rate from nucleoplasm into cytoplasm
#         kD: (cytoplasmic) RNA degradation rate
#         kch: RNA degradation rate at chromatin, set to zero
#         knp: RNA degradation rate in nucleoplasm, set to zero
#         km: mature RNA degradation rate, set to kD
#         kp: polysome RNA degradation rate, set to kD
#         Note: because km = kD = kp, kL drops out of equation (N7p193), so 
#         use cyto fraction for fitting kD. 
#         For fitting kL, use poly instead.
#     """
#     if k_para == None:
#         if len(k) == 3:#for Bayes Inference only relevant parameters are included, so not kL
#             kR = k[0]
#             kE = k[1] 
#             kD = k[2]    
#             kL = 1 #placeholder value, kL is not used as described above
#         else:#maintain backwards compatibility
#             kR = k[0]
#             kE = k[1]
#             kL = k[2] 
#             kD = k[3]
#     else:
#         kR = k_para[0]
#         kE = k_para[1]
#         kL = k_para[2]
#         kD = k[0]
#     kch = 0
#     knp = 0
#     km = kD
#     kp = kD
#     #Expressions get too large, so simplify
#     a = kch + kR
#     b = knp + kE
#     c = km + kL
#     d = kp + kL #NB: c == d when km = kD = kp
#     prefactor = b * c * kp / d 
#     Psi = ((a - b) * (c - a) * (kp - a))**(-1)
#     Phi = a / ((b - c) * b * (a - b) * (kp - b))
#     Theta =  a / (c * (b - c) * (c - a) * (kp - c))
#     Omega = (b * c * kp)**(-1)
#     Lambda = 1 + prefactor * ( (d - a) * Psi * np.exp(-a * t) + \
#              (d - b) * Phi * np.exp(-b * t) + \
#              (d - c) * Theta * np.exp(-c * t) \
#              - kL * (Psi + Phi + Theta + Omega) * np.exp(-kp * t) )
#     return Lambda

def lam_poly(k, t, k_para=None):
    """ Fit from nuc, cyto, poly compartments
        t: time
        k: rates list [kE, kL, kD]
        kE: RNA export rate from nucleus into cytoplasm
        kL: transition rate from cytoplasmic untranslated fraction into polysome fraction
        kD: (cytoplasmic) RNA degradation rate
    """
    if k_para == None:
        kE = k[0]
        kL = k[1]
        kD = k[2]
    else:
        kE = k_para[0]
        kL = k[0]
        kD = k_para[1]
    Phi = kD * (kD + kL) / ((kD + kL - kE) * (kD - kE))
    Psi = kE * kD / (kL * (kD + kL - kE))
    Omega = -kE * (kD + kL) / (kL * (kD - kE))
    Lambda = Phi * (1 - np.exp(-kE * t)) + \
             Psi * (1 - np.exp(-(kD + kL) * t)) + \
             Omega * (1 - np.exp(-kD * t))
    return Lambda

# def lam_poly_from_chr(k, t, k_para=None):
#     """ Fit from chr, nucpl, cyto, poly compartments
#         t: time
#         k: rates list [kR, kE, kL, kD]
#         kR: RNA Release rate from chromatin into nucleoplasm
#         kE: RNA export rate from nucleoplasm into cytoplasm
#         kL: transition rate from cytoplasmic untranslated fraction into polysome fraction
#         kD: (cytoplasmic) RNA degradation rate
#         kch: RNA degradation rate at chromatin, set to zero
#         knp: RNA degradation rate in nucleoplasm, set to zero
#         km: mature RNA degradation rate, set to kD.
#         kp: polysome RNA degradation rate, set to kD
#     """
#     if k_para == None:
#         kR = k[0]
#         kE = k[1]
#         kL = k[2] 
#         kD = k[3]
#     else:
#         kR = k_para[0]
#         kE = k_para[1]
#         kL = k[0]
#         kD = k_para[2]
#     kch = 0
#     knp = 0
#     km = kD
#     kp = kD
#     #Theta expression is too large, so simplify
#     a = kch + kR
#     b = knp + kE
#     c = km + kL
#     prefactor = b * c * kp 
#     Psi = ((a - b) * (c - a) * (kp - a))**(-1)
#     Phi = a / ((b - c) * b * (a - b) * (kp - b))
#     Theta =  a / (c * (b - c) * (c - a) * (kp - c))
#     Omega = prefactor**(-1)
#     Lambda = 1 + prefactor * ( Psi * np.exp(-a * t) + \
#              Phi * np.exp(-b * t) + \
#              Theta * np.exp(-c * t) \
#              - (Psi + Phi + Theta + Omega) * np.exp(-kp * t) )
#     return Lambda

# def lam_total(k, t, k_para=None):
#     """ Prediction from nuc, cyto, poly compartments
#         t: time
#         k: rates list [kE, kL, kD]
#         kE: RNA export rate from nucleus into cytoplasm
#         kD: (cytoplasmic) degradation rate
#     """
#     if k_para == None:
#         if len(k) == 2:#for Bayes Inference only relevant parameters are included, so not kL
#             kE = k[0] 
#             kD = k[1]          
#         else:#maintain backwards compatibility
#             kE = k[0]
#             kD = k[2]
#     else:
#         kE = k_para
#         kD = k[0]
#     if kD == kE:
#         Lambda = 1 - np.exp(-kD * t)
#     else: 
#         Phi = kD**2 / ((kD**2) - (kE**2))
#         Omega = -(kE**2) / ((kD**2) - (kE**2)) #1-Phi
#         Lambda = Phi * (1 - np.exp(-kE * t)) + \
#                  Omega * (1 - np.exp(-kD * t))
#     return Lambda

def lam_total_one_step(k, t, k_para=None):
    """ Fit total RNA degradation rate
        through a one step process
        t: time
        k: rates list [kD]
        kD: (total) degradation rate
    """
    kD = k[0]
    Lambda = 1 - np.exp(-kD * t)
    return Lambda

# def lam_total_from_chr(k, t, k_para=None):
#     """ Prediction from chr, nucpl, nuc, cyto, poly compartments
#         t: time
#         k: rates list [kR, kE, kL, kD]
#         kR: RNA Release rate from chromatin into nucleoplasm
#         kE: RNA export rate from nucleoplasm into cytoplasm
#         kL: transition rate from cytoplasmic untranslated fraction into polysome fraction
#         kD: (cytoplasmic) RNA degradation rate
#         kch: RNA degradation rate at chromatin, set to zero
#         knp: RNA degradation rate in nucleoplasm, set to zero
#         km: mature RNA degradation rate, set to kD
#         kp: polysome RNA degradation rate, set to kD
#         Note: because km = kD = kp, kL drops out of equation (N7p193).
#     """
#     if k_para == None:
#         kR = k[0]
#         kE = k[1]
#         kL = k[2] 
#         kD = k[3]
#     else:
#         kR = k_para[0]
#         kE = k_para[1]
#         kL = k_para[2]
#         kD = k_para[3]
#     kch = 0
#     knp = 0
#     km = kD
#     kp = kD
#     #Expressions get too large, so simplify
#     a = kch + kR
#     b = knp + kE
#     c = km + kL
#     d = kp + kL #NB: c == d when km = kD = kp
#     prefactor = b * c * kp / (b * c * kp + c * kR * kp + kR * kE * d) 
#     Psi = ((a - b) * (c - a) * (kp - a))**(-1)
#     Phi = a / ((b - c) * b * (a - b) * (kp - b))
#     Theta =  a / (c * (b - c) * (c - a) * (kp - c))
#     Omega = (b * c * kp)**(-1)
#     Lambda = 1 + prefactor * ( \
#         ((kR + b - a) / (a - b) + kR * kE * (d - a) * Psi) * np.exp(-a * t) + \
#         (-kR * a / (b * (a - b)) + kR * kE * (d - b) * Phi) * np.exp(-b * t) + \
#         kR * kE * (d - c) * Theta * np.exp(-c * t) \
#         - kR * kE * kL * (Psi + Phi + Theta + Omega) * np.exp(-kp * t) )
#     return Lambda

def one_error_background(params, TCdist, fixed_params=None):
    pe1 = params['Error_1']    
    (N,K) = TCdist.shape    
    model1 = np.zeros((N,K))
    
    for n in range(N):
        TCsum = sum(TCdist[n,])
        for k in range(K):
            if n >= k:
                model1[n, k] = binom.pmf(k,n,pe1) * TCsum 
    return model1

def two_error_background(params, TCdist, fixed_params=None):
    pe1 = params['Error_1']
    pe2 = params['Error_2']
    pi = params['frac_err']
    
    (N,K) = TCdist.shape    
    model2 = np.zeros((N,K))
    
    for n in range(N):
        TCsum = sum(TCdist[n,])
        for k in range(K):
            if n >= k:
                model2[n, k] = (pi * binom.pmf(k,n,pe1) + (1-pi) * binom.pmf(k,n,pe2) ) * TCsum 
    return model2

def TC_rate_from_background(params, TCdist, fixed_params):
    pe1 = fixed_params['Error_1']
    pe2 = fixed_params['Error_2']
    pi_err = fixed_params['frac_err']
    pc = params['TC']
    pi = params['frac']
    
    (N,K) = TCdist.shape    
    model_TC = np.zeros((N,K))
    
    for n in range(N):
        TCsum = sum(TCdist[n,])
        for k in range(K):
            if n >= k:
                model_TC[n, k] = ( pi * binom.pmf(k,n,pc) + 
                 (1 - pi) * (pi_err * binom.pmf(k,n,pe1) + (1 - pi_err) * binom.pmf(k,n,pe2)) ) * TCsum 
    return model_TC

def TCconv(k, t, k_para=None):
    """ 4sU dynamics
        t: time
        kin: 4sU in rate
        kout: 4sU out rate derivation N6p184: TC = 1/t * int_{0}^{t} L(s)ds
    """
    if type(k) == lmfit.parameter.Parameters:
        kin = k['kin']
        kout = k['kout']
    else:
        kin = k[0]
        kout = k[1]
    TC = kin / kout *( 1 + (np.exp(-kout * t) - 1) / (kout * t) )
    return TC

def pc_instant(k, t, k_para=None):
    """ 4sU dynamics
        t: time
        kin: 4sU in rate
        kout: 4sU out rate derivation N6p184: TC = int_{0}^{t} L(s)ds
    """
    if type(k) == lmfit.parameter.Parameters:
        kin = k['kin']
        kout = k['kout']
    else:
        kin = k[0]
        kout = k[1]
    pc = kin / kout * ( 1 - np.exp(-kout * t) )
    return pc