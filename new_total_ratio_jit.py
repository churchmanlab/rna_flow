#!/usr/bin/env python
# coding: utf-8

# ## new_total_ratio.py
# Author: Robert Ietswaart
# Date: 20210708
# License: BSD2.  
# Python v3.7.4 
# 
# Subcellular Timelapse Seq: fit time scales to his subcellular timelapse seq using observed new to total ratio (Lambda)
# of the various subcellular RNA fractions
# See Notebook N7p51-56,p174-179,p181-191, N8p96-109
# Comparison with new_total_ratio.py:
# -jit additions
# -lam_nuc_from_chr (with nucdeg): N8p201 and N7 as above.
# -lam_total used as lam_nuc_from_chr (no nucdeg): N8p193-6 for derivation.
# -lam_cyto_from_chr is modified to explicitly use the closed form that is independent of kL, derivs: N8p201, N7p193 and N8p-3 
# -lam_total_from_chr idem: N7p193 and N8p-3, N8p14-15, and N8p201 for form with kND 
# -lam_poly_from_chr is modified to explicitly use the simplified closed form, derivations: N8p1, N8p21-22
# -lam_poly_from_nuc idem: N8p18
# -Nuclear degradation model: N8p96-108
# -limit cases where rates differ less than eps: N8p8, N8p18-20, N8p109, N8p205

import numpy as np
from numba import jit, float64

@jit(float64(float64, float64), nopython=True, cache=True)
def lam_1var(a, t):
    """ Derivative of Lambda with respect to only var a at t
        t: time point (float)  
        a: rate a in domain [0,inf)
        der: derivative, function on a domain
    """
    Lambda = 1 - np.exp(-a * t)
    if Lambda < 0:#numerical under/overflow
        Lambda = 0
    return Lambda

@jit(float64(float64, float64), nopython=True, cache=True)
def lam_chr(a, t):
    """ Fit based on chr compartment
        t: time
        kR: RNA Release rate from chromatin into nucleoplasm
        knd_ch: RNA degradation rate at chromatin
        a: kR + knd_ch, the chromatin residence rate
    """
    return lam_1var(a, t)

@jit(float64(float64, float64), nopython=True, cache=True)
def lam_nuc(a, t):
    """ Fit based on only nuc compartment data.
        t: time
        a: nuclear residence rate
    """
    return lam_1var(a, t)

@jit(float64(float64, float64, float64, float64), nopython=True, cache=True)
def lam_nuc_from_chr(a, b, kND, t):
    """ Prediction as sum of chr and nucpl compartments
        t: time
        kR: RNA Release rate from chromatin into nucleoplasm
        kE: RNA Export rate from nucleoplasm into cytoplasm
        kND: Nuclear RNA Degradation rate (assumed equal on chromatin and nucleoplasm)
        a: kR + kND 
        b: kE + kND (so a,b >= kND)
    """
    eps = 1e-16
    if np.absolute(((a + b - kND) * (a - b))) <= eps:
        Lambda = 1 - (1 + b * (b - kND) * t / (2 * b - kND)) * np.exp(-b * t)
#         print('warning: nuc_from_chr exact limit case, lam:', Lambda, a, b, kND, t)
    else:
        prefactor = ((a + b - kND) * (a - b))**(-1)
        Phi = b * (b - kND)
        Omega = -a * (a - kND)
        Lambda = 1 + prefactor * (Phi * np.exp(-a * t) + Omega * np.exp(-b * t))     
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_nuc_from_chr<0: ', Lambda, a, b, kND, t)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64, float64), nopython=True, cache=True)
def lam_cyto(a, kD, t):
    """ Fit based on nuc, cyto compartment data.
        t: time
        a: nuclear residence rate
        kD: cytoplasmic RNA degradation rate
    """
    eps = 1e-16
    if np.absolute(kD - a) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t)
#         print('cyto exact limit case, lam:', Lambda, a, kD, t)
    else:
        Phi = kD / (kD - a)
        Omega = -a / (kD - a) #1-Phi
        Lambda = Phi * (1 - np.exp(-a * t)) + \
                 Omega * (1 - np.exp(-kD * t))      
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_cyto<0: ', Lambda, a, kD, t)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64, float64, float64), nopython=True, cache=True)
def lam_cyto_from_chr(a, b, kD, t):
    """ Fit based on chr, nuc and cyto compartments
        t: time
        kR: RNA Release rate from chromatin into nucleoplasm
        kE: RNA export rate from nucleoplasm into cytoplasm
        kD: cytoplasmic RNA Degradation rate
        knd_ch: Nuclear RNA Degradation rate on chromatin
        knd_np: Nuclear RNA Degradation rate in nucleoplasm
        a: kR + knd_ch
        b: kE + knd_np
    """
    eps = 1e-16
    if np.absolute(a - b) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t)
#         print('warning: cyto_from_chr approx limit case1, lam:', Lambda, a, b, kD, t)        
    elif np.absolute(kD - a) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t)
#         print('warning: cyto_from_chr approx limit case2, lam:', Lambda, a, b, kD, t)
    elif np.absolute(kD - b) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t)
#         print('warning: cyto_from_chr approx limit case3, lam:', Lambda, a, b, kD, t)
    else:
        Psi = b * kD * ((a - b) * (kD - a))**(-1)
        Phi = -a * kD * ((a - b) * (kD - b))**(-1)
        Omega = -a * b * ((kD - a) * (kD - b))**(-1)
        Lambda = 1 + Psi * np.exp(-a * t) + \
                 Phi * np.exp(-b * t) + \
                 Omega * np.exp(-kD * t)    
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_cyto_from_chr<0: ', Lambda, kR, kE, kD, t)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64, float64, float64), nopython=True, cache=True)
def lam_poly(a, kD, kL, t):
    """ Fit from nuc, cyto, poly compartments
        t: time
        a: Nuclear residence rate
        kD: Cytoplasmic RNA degradation rate
        kL: transition rate from cytoplasmic untranslated fraction into polysome fraction
    """
    eps = 1e-12 #takes 40 mins to evaluate, with 1e-16: lambda turns negative, 1e-6 no improvement
    c = kD + kL
    if np.absolute(kD - a) <= eps: 
        Lambda = 1 - (1 + a**2 / (kL * (c - a)) - a * c * t / kL) * np.exp(-a * t) + \
            -(a**2) / (kL * (c - a)) * np.exp(-c * t)
#         if Lambda < 0: 
#             print('poly exact limit case1, lam:', Lambda, a, kD, kL, c, t)
    elif np.absolute(c - a) <= eps: 
        Lambda = 1 - (1 - a * kD * t / kL + a**2 / (kL * (kD - a))) * np.exp(-a * t) + \
            a**2 / (kL * (kD - a)) * np.exp(-kD * t)
#         if Lambda < 0: 
#             print('poly exact limit case2, lam:', Lambda, a, kD, kL, c, t)
    else:
        Psi = -kD * c / ((c - a) * (kD - a))
        Phi = -a * kD / (kL * (c - a))
        Omega = a * c / (kL * (kD - a))
        Lambda = 1 + Psi * np.exp(-a * t) + \
                 Phi * np.exp(-c * t) + \
                 Omega * np.exp(-kD * t)      
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_poly<0: ', Lambda, a, kD, kL, c, t)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64), nopython=True, cache=True)
def lam_total_one_step(kD, t):
    """ Fit total RNA degradation rate
        through a one step process
        t: time
        kD: total RNA Degradation rate, aka whole cell turnover rate
    """
    return lam_1var(kD, t)

@jit(float64(float64, float64, float64), nopython=True, cache=True)
def lam_total(a, kD, t):
    """ Prediction from nuc, cyto compartments
        t: time
        a: nuclear residence rate
        kD: cytoplasmic Degradation rate
    """
    eps = 1e-16
    if np.absolute(kD - a) <= eps:
        Lambda = 1 - (1 + kD * t / 2) * np.exp(-kD * t)
#         print('warning: total exact limit case, lam:', Lambda, a, t)
    else: 
        Phi = kD**2 / ((kD**2) - (a**2))
        Omega = -(a**2) / ((kD**2) - (a**2)) #1-Phi
        Lambda = Phi * (1 - np.exp(-a * t)) + \
                 Omega * (1 - np.exp(-kD * t))   
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_total<0: ', Lambda, a, kD, t)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64, float64, float64, float64), nopython=True, cache=True)
def lam_total_from_chr(a, b, kD, kND, t):
    """ Prediction from chr, nucpl, nuc, cyto, poly compartments
        t: time
        kR: RNA Release rate from chromatin into nucleoplasm
        kE: RNA Export rate from nucleoplasm into cytoplasm
        kD: cytoplasmic RNA Degradation rate
        kND: Nuclear RNA Degradation rate (assumed equal on chromatin and nucleoplasm)
        a: kR + kND
        b: kE + kND
    """  
    eps = 1e-16 #1e-12 does not improve runtime with limit approximations
    if np.absolute(a - b) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t)
#         print('warning: total_from_chr approx limit case1, lam:', Lambda, a, b, kD, kND, t)        
    elif np.absolute(kD - a) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t)
#         print('warning: total_from_chr approx limit case2, lam:', Lambda, a, b, kD, kND, t)
    elif np.absolute(kD - b) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t)
#         print('warning: total_from_chr approx limit case3, lam:', Lambda, a, b, kD, kND, t)
    else:
        prefactor = b * kD * (kD * (a + b - kND) + (a - kND) * (b - kND))**(-1)
        Psi = (b - kND) * (kD - kND) * ((a - b) * (kD - a))**(-1)
        Phi = - a * (a - kND) * (kD - kND) * (b * (a - b) * (kD - b))**(-1)
        Omega = -a * (a - kND) * (b - kND) * (kD * (kD - a) * (kD - b))**(-1)
        Lambda = 1 + prefactor * ( \
                 Psi * np.exp(-a * t) + \
                 Phi * np.exp(-b * t) + \
                 Omega * np.exp(-kD * t) )     
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_total_from_chr<0: ', Lambda, a, b, kD, kND, t)
        Lambda = 0
    return Lambda

#3 compartment Nuclear Degradation model
@jit(float64(float64, float64, float64, float64), nopython=True, cache=True)
def lam_total_from_nucdeg(a, kD, k_, t):
    """ Fit total RNA
        t: time point (float) 
        kND: nuclear degradation rate
        k_: nuclear exit flow rate (not to be used further, 
        not to be interpreted as nuclear export rate)
        a: kND + k_: nuclear residence rate (kN)
        kD: cytoplasmic residence (degradation) rate
    """
    eps = 1e-16
    if np.absolute(a - kD) <= eps:
        Lambda = 1 - (1 + k_ * a * t * (k_ + kD)**(-1)) * np.exp(-a * t)
#         print('warning: total_nucdeg exact limit case, Lambda, a, kD, k_, t)
    else: 
        Psi = ((a - kD) * (kD + k_))**(-1)
        Lambda = 1 - Psi * kD * (a - kD - k_) * np.exp(-a * t) \
                 - Psi * k_ * a * np.exp(-kD * t)
#         print('lam: total_nucdeg', Lambda, a, kD, k_, t)
    if Lambda < 0:#numerical under/overflow
#         print('warning: total_nucdeg, lam<0', Lambda, a, kD, k_, t)
        Lambda = 0
    return Lambda