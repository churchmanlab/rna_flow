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
# -lam_cyto_from_chr is modified to explicitly use the closed form that is independent of kL, derivations: N7p193 and N8p-3
# -lam_total_from_chr idem: N7p193 and N8p-3, N8p14-15
# -lam_poly_from_chr is modified to explicitly use the simplified closed form, derivations: N8p1, N8p21-22
# -lam_poly_from_nuc idem: N8p18
# -Nuclear degradation model: N8p96-108
# -limit cases where rates differ less than eps: N8p8, N8p18-20, N8p109

import numpy as np
from numba import jit, float64

@jit(float64(float64, float64), nopython=True, cache=True)
def lam_chr(kR, t):
    """ Fit based on chr compartment
        t: time
        kR: RNA Release rate from chromatin into nucleoplasm
        kch: RNA degradation rate at chromatin, set to zero
    """
    kch = 0
    Lambda = 1 - np.exp(-(kR + kch) * t)
    if Lambda < 0:#numerical under/overflow
#         print('lam_chr<0: ', Lambda, kR, t)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64), nopython=True, cache=True)
def lam_nuc(a, t):
    """ Fit based on only nuc compartment data.
        t: time
        a: Nuclear residence rate
    """
    Lambda = 1 - np.exp(-a * t)
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_nuc<0: ', Lambda, a)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64, float64), nopython=True, cache=True)
def lam_nuc_from_chr(kR, kE, t):
    """ Prediction as sum of chr and nucpl compartments
        t: time
        kR: RNA Release rate from chromatin into nucleoplasm
        kE: RNA export rate from nucleoplasm into cytoplasm
        kch: RNA degradation rate at chromatin, set to zero
        knp: RNA degradation rate in nucleoplasm, set to zero
    """
    eps = 1e-16
    kch = 0
    knp = 0
    a = kch + kR
    b = knp + kE
    if np.absolute(a - b) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t)#limit approximation, but nth order effect
#         print('warning: nuc_from_chr limit case, lam:', Lambda, kR, kE, t)
    else:
        prefactor = ((kR + b) * (a - b))**(-1)
        Phi = -b * (kR + b - a)
        Omega = kR * a
        Lambda = prefactor * ( Phi * (1 - np.exp(-a * t)) + \
                 Omega * (1 - np.exp(-b * t)) )     
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_nuc_from_chr<0: ', Lambda, kR, kE, t)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64, float64), nopython=True, cache=True)
def lam_cyto(a, kD, t):
    """ Fit based on nuc, cyto compartment data.
        t: time
        a: Nuclear residence rate
        kD: (cytoplasmic) RNA degradation rate
    """
    eps = 1e-16
    if np.absolute(kD - a) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t)
#         print('cyto limit case: exact, lam:', Lambda, a, kD, t)
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
def lam_cyto_from_chr(kR, kE, kD, t):
    """ Fit based on chr, nucpl and cyto compartments
        t: time
        kR: RNA Release rate from chromatin into nucleoplasm
        kE: RNA export rate from nucleoplasm into cytoplasm
        kD: (cytoplasmic) RNA degradation rate
        kch: RNA degradation rate at chromatin, set to zero
        knp: RNA degradation rate in nucleoplasm, set to zero
        km: mature RNA degradation rate, set to kD
        kp: polysome RNA degradation rate, set to kD
        Note: because km = kD = kp, kL drops out of equation (N7p193), so 
        use cyto fraction for fitting kD. 
        For fitting kL, use poly instead.
    """
    eps = 1e-16
    kch = 0
    knp = 0
    a = kch + kR
    b = knp + kE
    if np.absolute(a - b) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation, but nth order effect
#         print('warning: cyto_from_chr limit case1, lam:', Lambda, kR, kE, kD, t)        
    elif np.absolute(kD - a) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation, but nth order effect
#         print('warning: cyto_from_chr limit case2, lam:', Lambda, kR, kE, kD, t)
    elif np.absolute(kD - b) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation, but nth order effect
#         print('warning: cyto_from_chr limit case3, lam:', Lambda, kR, kE, kD, t)
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
def lam_poly(kE, kD, kL, t):
    """ Fit from nuc, cyto, poly compartments
        t: time
        kE: RNA export rate from nucleus into cytoplasm
        kD: (cytoplasmic) RNA degradation rate
        kL: transition rate from cytoplasmic untranslated fraction into polysome fraction
    """
    eps = 1e-12 #takes 40 mins to evaluate, with 1e-16: lambda turns negative, 1e-6 no improvement
    c = kD + kL
    if np.absolute(kD - kE) <= eps: 
        Lambda = 1 - (1 + kE**2 / (kL * (c - kE)) - kE * c * t / kL) * np.exp(-kE * t) + \
            -(kE**2) / (kL * (c - kE)) * np.exp(-c * t)
#         if Lambda < 0: 
#             print('poly limit case1, lam: exact')#, Lambda, kE, kD, kL, c, t)
    elif np.absolute(c - kE) <= eps: 
        Lambda = 1 - (1 - kE * kD * t / kL + kE**2 / (kL * (kD - kE))) * np.exp(-kE * t) + \
            kE**2 / (kL * (kD - kE)) * np.exp(-kD * t)
#         if Lambda < 0: 
#             print('poly limit case2, lam: exact')#, Lambda, kE, kD, kL, c, t)
    else:
        Psi = -kD * c / ((c - kE) * (kD - kE))
        Phi = -kE * kD / (kL * (c - kE))
        Omega = kE * c / (kL * (kD - kE))
        Lambda = 1 + Psi * np.exp(-kE * t) + \
                 Phi * np.exp(-c * t) + \
                 Omega * np.exp(-kD * t)      
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_poly<0: ', Lambda, kE, kD, kL, c, t)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64, float64, float64, float64), nopython=True, cache=True)
def lam_poly_from_chr(kR, kE, kD, kL, t):
    """ Fit from chr, nucpl, cyto, poly compartments
        t: time
        kR: RNA Release rate from chromatin into nucleoplasm
        kE: RNA export rate from nucleoplasm into cytoplasm
        kD: (cytoplasmic) RNA degradation rate
        kL: transition rate from cytoplasmic untranslated fraction into polysome fraction
        kch: RNA degradation rate at chromatin, set to zero
        knp: RNA degradation rate in nucleoplasm, set to zero
        km: mature RNA degradation rate, set to kD.
        kp: polysome RNA degradation rate, set to kD
    """
    eps = 1e-12 #1e-6 causes case2 to occur
    kch = 0
    knp = 0
    a = kch + kR
    b = knp + kE
    c = kD + kL
    if np.absolute(a - b) <= eps: 
        Theta = -(b**2) * kD * ((b - c)**2 * (kD - c))**(-1)
        Xi = a * b**2 * ((kD - c) * (kD - b)**2)**(-1) #-(1+Psi+Phi+Theta)
        Phi = 1 + b * c * kD * t * ((b - c) * (kD - b))**(-1) + Theta + Xi
        Lambda = 1 - Phi * np.exp(-b * t) + Theta * np.exp(-c * t) + Xi * np.exp(-kD * t)
#         if Lambda < 0:
#             print('poly_from_chr limit case1, exact lam: ')#, Lambda, kR, kE, kD, kL, c, t)        
    elif np.absolute(a - c) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation
#         if Lambda < 0:
#         print('warning: poly_from_chr limit case2, lam:')#, Lambda, kR, kE, kD, kL, c, t)
    elif np.absolute(b - c) <= eps:
        Psi = -(b**2) * kD * ((a - b)**2 * (kD - a))**(-1)
        Xi = a * b**2 * ((kD - a) * (kD - b)**2)**(-1) #-(1+Psi+Phi+Theta)
        Phi = 1 + Psi + a * b * kD * t * ((a - b) * (kD - b))**(-1) + Xi
        Lambda = 1 + Psi * np.exp(-a * t) - Phi * np.exp(-b * t) + Xi * np.exp(-kD * t)
#         if Lambda < 0:
#         print('poly_from_chr limit case3, exact lam:')#, Lambda, kR, kE, kD, kL, c, t)
    elif np.absolute(kD - a) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation
#         if Lambda < 0:
#         print('warning: poly_from_chr limit case4, lam:')#, Lambda, kR, kE, kD, kL, c, t)
    elif np.absolute(kD - b) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation
#         if Lambda < 0:
#         print('warning: poly_from_chr limit case5, lam:')#, Lambda, kR, kE, kD, kL, c, t)
    elif np.absolute(kD - c) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation
#         if Lambda < 0:
#         print('warning: poly_from_chr limit case6, lam:')#, Lambda, kR, kE, kD, kL, c, t)
    else:    
        Psi = b * c * kD * ((a - b) * (c - a) * (kD - a))**(-1)
        Phi = a * c * kD * ((a - b) * (b - c) * (kD - b))**(-1)
        Theta = a * b * kD * ((b - c) * (c - a) * (kD - c))**(-1)
#         Xi = a * b * c * ((kD - a) * (kD - b) * (kD - c))**(-1) #-(1+Psi+Phi+Theta)
        Lambda = 1 + Psi * np.exp(-a * t) + \
                 Phi * np.exp(-b * t) + \
                 Theta * np.exp(-c * t) + \
                 -(1 + Psi + Phi + Theta) * np.exp(-kD * t) #Xi * np.exp(-kD * t)
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_poly_from_chr<0: ', Lambda, kR, kE, kD, kL, c, t)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64, float64), nopython=True, cache=True)
def lam_total(a, kD, t):
    """ Prediction from nuc, cyto, poly compartments
        t: time
        a: Nuclear residence rate
        kD: (cytoplasmic) degradation rate
    """
    eps = 1e-16
    if np.absolute(kD - a) <= eps:
        Lambda = 1 - np.exp(-a * t) - a * t * np.exp(-a * t) #limit approximation, but nth order effect
#         print('warning: total limit case, lam:', Lambda, kD, kE, t)
    else: 
        Phi = kD**2 / ((kD**2) - (a**2))
        Omega = -(a**2) / ((kD**2) - (a**2)) #1-Phi
        Lambda = Phi * (1 - np.exp(-a * t)) + \
                 Omega * (1 - np.exp(-kD * t))   
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_total<0: ', Lambda, a, kD, t)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64), nopython=True, cache=True)
def lam_total_one_step(kD, t):
    """ Fit total RNA degradation rate
        through a one step process
        t: time
        kD: (total) degradation rate
    """
    Lambda = 1 - np.exp(-kD * t)
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_total_one_step<0: ', Lambda, kD)
        Lambda = 0
    return Lambda

@jit(float64(float64, float64, float64, float64), nopython=True, cache=True)
def lam_total_from_chr(kR, kE, kD, t):
    """ Prediction from chr, nucpl, nuc, cyto, poly compartments
        t: time
        kR: RNA Release rate from chromatin into nucleoplasm
        kE: RNA export rate from nucleoplasm into cytoplasm
        kD: (cytoplasmic) RNA degradation rate
        kch: RNA degradation rate at chromatin, set to zero
        knp: RNA degradation rate in nucleoplasm, set to zero
        Note: because km = kD = kp, kL drops out of equation (N7p193), 
        so kL is not an argument.
    """  
    eps = 1e-16
    kch = 0
    knp = 0
    a = kch + kR
    b = knp + kE
    if np.absolute(a - b) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation, but nth order effect
#         print('warning: total_from_chr limit case1, lam:', Lambda, kR, kE, kD, t)        
    elif np.absolute(kD - a) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation, but nth order effect
#         print('warning: total_from_chr limit case2, lam:', Lambda, kR, kE, kD, t)
    elif np.absolute(kD - b) <= eps:
        Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation, but nth order effect
#         print('warning: total_from_chr limit case3, lam:', Lambda, kR, kE, kD, t)
    else:
        prefactor = b * kD * (b * kD + kR * kD + kR * kE)**(-1)
        Psi = (kR + b - a) * (a - b)**(-1) + kR * kE * ((a - b) * (kD - a))**(-1)
        Phi = -kR * a * (b * (a - b))**(-1) - a * kR * kE * (b * (a - b) * (kD - b))**(-1)
        Omega = -a * kR * kE * (kD * (kD - a) * (kD - b))**(-1)
        Lambda = 1 + prefactor * ( \
                 Psi * np.exp(-a * t) + \
                 Phi * np.exp(-b * t) + \
                 Omega * np.exp(-kD * t) )     
    if Lambda < 0:#numerical under/overflow
#         print('warning: lam_total_from_chr<0: ', Lambda, kR, kE, kD, t)
        Lambda = 0
    return Lambda

#Nuclear Degradation model
@jit(float64(float64, float64, float64, float64), nopython=True, cache=True)
def lam_total_from_nucdeg(a, kD, kE, t):
    """ Fit total RNA
        t: time point (float) 
        kN: nuclear degradation rate
        kE: nuclear export rate
        a: kN + kE: nuclear residence rate
        kD: cytoplasmic residence (degradation) rate
    """
    eps = 1e-16
    if np.absolute(a - kD) <= eps:
        Lambda = 1 - (1 + kE * a * t * (kE + kD)**(-1)) * np.exp(-a * t)#exact limit
#         print('warning: total_nucdeg limit case, Lambda, a, kD, kE, t)
    else: 
        Psi = ((a - kD) * (kD + kE))**(-1)
        Lambda = 1 - Psi * kD * (a - kD - kE) * np.exp(-a * t) \
                 - Psi * kE * a * np.exp(-kD * t)
#         print('lam: total_nucdeg', Lambda, a, kD, kE, t)
    if Lambda < 0:#numerical under/overflow
#         print('warning: total_nucdeg, lam<0', Lambda, a, kD, kE, t)
        Lambda = 0
    return Lambda




