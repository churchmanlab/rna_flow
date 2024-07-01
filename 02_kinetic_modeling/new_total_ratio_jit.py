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


# #OLD: uncomment for backwards compatibility with older versions of scripts
# @jit(float64(float64, float64), nopython=True, cache=True)
# def lam_nucdeg(a, t):
#     """ Fit nuclear fraction
#         t: time point (float)
#         kN: Nuclear degradation rate
#         kE: Nuclear export rate
#         a: kN + kE
#     """
#     return lam_total_one_step(a, t)

# @jit(float64(float64, float64, float64), nopython=True, cache=True)
# def lam_cyto_from_nucdeg(a, kD, t):
#     """ Fit cytoplasmic fraction
#         t: time point (float) 
#         kN: nuclear degradation rate
#         kE: nuclear export rate
#         a: kN + kE
#         kD: cytoplasmic residence (degradation) rate
#     """
#     return lam_cyto(a, kD, t)


# Not used
# @jit(float64(float64, float64, float64), nopython=True, cache=True)
# def lam_nucpl(kR, kE, t):
#     """ Fit based on chr and nucpl compartments
#         t: time
#         kR: RNA Release rate from chromatin into nucleoplasm
#         kE: RNA export rate from nucleoplasm into cytoplasm
#         kch: RNA degradation rate at chromatin, set to zero
#         knp: RNA degradation rate in nucleoplasm, set to zero
#     """
#     eps = 1e-16
#     kch = 0
#     knp = 0
#     a = kch + kR
#     b = knp + kE
#     if np.absolute(a - b) <= eps:
#         Lambda = 1 - (1 + a * t) * np.exp(-a * t)
# #         print('nucpl limit case: exact, lam:', Lambda, kR, kE, t)
#     else:
#         Phi = -a / (a - b)
#         Omega = b / (a - b) #1-Phi
#         Lambda = 1 + Phi * np.exp(-b * t) + \
#                  Omega * np.exp(-a * t)
#     if Lambda < 0:#numerical under/overflow
# #         print('warning: lam_nucpl<0: ', Lambda, kR, kE, t)
#         Lambda = 0
#     return Lambda



# @jit(float64(float64, float64, float64, float64, float64), nopython=True, cache=True)
# def lam_poly_from_chr(kR, kE, kD, kL, t):
#     """ Fit from chr, nucpl, cyto, poly compartments
#         t: time
#         kR: RNA Release rate from chromatin into nucleoplasm
#         kE: RNA export rate from nucleoplasm into cytoplasm
#         kD: (cytoplasmic) RNA degradation rate
#         kL: transition rate from cytoplasmic untranslated fraction into polysome fraction
#         kch: RNA degradation rate at chromatin, set to zero
#         knp: RNA degradation rate in nucleoplasm, set to zero
#         km: mature RNA degradation rate, set to kD.
#         kp: polysome RNA degradation rate, set to kD
#     """
#     eps = 1e-12 #1e-6 causes case2 to occur
#     kch = 0
#     knp = 0
#     a = kch + kR
#     b = knp + kE
#     c = kD + kL
#     if np.absolute(a - b) <= eps: 
#         Theta = -(b**2) * kD * ((b - c)**2 * (kD - c))**(-1)
#         Xi = a * b**2 * ((kD - c) * (kD - b)**2)**(-1) #-(1+Psi+Phi+Theta)
#         Phi = 1 + b * c * kD * t * ((b - c) * (kD - b))**(-1) + Theta + Xi
#         Lambda = 1 - Phi * np.exp(-b * t) + Theta * np.exp(-c * t) + Xi * np.exp(-kD * t)
# #         if Lambda < 0:
# #             print('poly_from_chr limit case1, exact lam: ')#, Lambda, kR, kE, kD, kL, c, t)        
#     elif np.absolute(a - c) <= eps:
#         Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation
# #         if Lambda < 0:
# #         print('warning: poly_from_chr limit case2, lam:')#, Lambda, kR, kE, kD, kL, c, t)
#     elif np.absolute(b - c) <= eps:
#         Psi = -(b**2) * kD * ((a - b)**2 * (kD - a))**(-1)
#         Xi = a * b**2 * ((kD - a) * (kD - b)**2)**(-1) #-(1+Psi+Phi+Theta)
#         Phi = 1 + Psi + a * b * kD * t * ((a - b) * (kD - b))**(-1) + Xi
#         Lambda = 1 + Psi * np.exp(-a * t) - Phi * np.exp(-b * t) + Xi * np.exp(-kD * t)
# #         if Lambda < 0:
# #         print('poly_from_chr limit case3, exact lam:')#, Lambda, kR, kE, kD, kL, c, t)
#     elif np.absolute(kD - a) <= eps:
#         Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation
# #         if Lambda < 0:
# #         print('warning: poly_from_chr limit case4, lam:')#, Lambda, kR, kE, kD, kL, c, t)
#     elif np.absolute(kD - b) <= eps:
#         Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation
# #         if Lambda < 0:
# #         print('warning: poly_from_chr limit case5, lam:')#, Lambda, kR, kE, kD, kL, c, t)
#     elif np.absolute(kD - c) <= eps:
#         Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation
# #         if Lambda < 0:
# #         print('warning: poly_from_chr limit case6, lam:')#, Lambda, kR, kE, kD, kL, c, t)
#     else:    
#         Psi = b * c * kD * ((a - b) * (c - a) * (kD - a))**(-1)
#         Phi = a * c * kD * ((a - b) * (b - c) * (kD - b))**(-1)
#         Theta = a * b * kD * ((b - c) * (c - a) * (kD - c))**(-1)
# #         Xi = a * b * c * ((kD - a) * (kD - b) * (kD - c))**(-1) #-(1+Psi+Phi+Theta)
#         Lambda = 1 + Psi * np.exp(-a * t) + \
#                  Phi * np.exp(-b * t) + \
#                  Theta * np.exp(-c * t) + \
#                  -(1 + Psi + Phi + Theta) * np.exp(-kD * t) #Xi * np.exp(-kD * t)
#     if Lambda < 0:#numerical under/overflow
# #         print('warning: lam_poly_from_chr<0: ', Lambda, kR, kE, kD, kL, c, t)
#         Lambda = 0
#     return Lambda

# @jit(float64(float64, float64, float64, float64), nopython=True, cache=True)
# def lam_total_from_chr(kR, kE, kD, t):
#     """ Prediction from chr, nucpl, nuc, cyto, poly compartments
#         t: time
#         kR: RNA Release rate from chromatin into nucleoplasm
#         kE: RNA export rate from nucleoplasm into cytoplasm
#         kD: Cytoplasmic RNA degradation rate
#     """  
#     eps = 1e-16
#     kch = 0
#     knp = 0
#     a = kch + kR
#     b = knp + kE
#     if np.absolute(a - b) <= eps:
#         Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation, but nth order effect
# #         print('warning: total_from_chr limit case1, lam:', Lambda, kR, kE, kD, t)        
#     elif np.absolute(kD - a) <= eps:
#         Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation, but nth order effect
# #         print('warning: total_from_chr limit case2, lam:', Lambda, kR, kE, kD, t)
#     elif np.absolute(kD - b) <= eps:
#         Lambda = 1 - (1 + a * t) * np.exp(-a * t) #limit approximation, but nth order effect
# #         print('warning: total_from_chr limit case3, lam:', Lambda, kR, kE, kD, t)
#     else:
#         prefactor = b * kD * (b * kD + kR * kD + kR * kE)**(-1)
#         Psi = (kR + b - a) * (a - b)**(-1) + kR * kE * ((a - b) * (kD - a))**(-1)
#         Phi = -kR * a * (b * (a - b))**(-1) - a * kR * kE * (b * (a - b) * (kD - b))**(-1)
#         Omega = -a * kR * kE * (kD * (kD - a) * (kD - b))**(-1)
#         Lambda = 1 + prefactor * ( \
#                  Psi * np.exp(-a * t) + \
#                  Phi * np.exp(-b * t) + \
#                  Omega * np.exp(-kD * t) )     
#     if Lambda < 0:#numerical under/overflow
# #         print('warning: lam_total_from_chr<0: ', Lambda, kR, kE, kD, t)
#         Lambda = 0
#     return Lambda


#implementation of nuclear degradation model
# @jit(float64(float64, float64, float64), nopython=True, cache=True)
# def lam_nucdeg(kN, kE, t):
#     """ Fit nuclear fraction
#         t: time point (float)
#         kN: Nuclear degradation rate
#         kE: Nuclear export rate
#     """
#     return lam_total_one_step((kN + kE), t)

# @jit(float64(float64, float64, float64, float64), nopython=True, cache=True)
# def lam_cyto_from_nucdeg(kN, kE, kD, t):
#     """ Fit cytoplasmic fraction
#         t: time point (float) 
#         kN: nuclear degradation rate
#         kE: nuclear export rate
#         kD: cytoplasmic residence (degradation) rate
#     """
#     return lam_cyto((kN + kE), kD, t)

# @jit(float64(float64, float64, float64, float64), nopython=True, cache=True)
# def lam_total_from_nucdeg(kN, kE, kD, t):
#     """ Fit total RNA
#         t: time point (float) 
#         kN: nuclear degradation rate
#         kE: nuclear export rate
#         kD: cytoplasmic residence (degradation) rate
#     """
#     eps = 1e-16
#     a = kN + kE
#     Ap = kD * (kN - kD)# = kD * (a - (kD + kE))
#     Bp = kE * a
#     Lambda = 1 - Ap / (Ap + Bp) * np.exp(-a * t) \
#              - Bp / (Ap + Bp) * np.exp(-kD * t)
#     if Lambda < 0:#numerical under/overflow
#         Lambda = 0
#     return Lambda
    

# #Mitoribo model: transferred to RNAdecay project scripts
# @jit(float64(float64, float64, float64), nopython=True, cache=True)
# def lam_ribo(kL, kD, t):
#     """ Fit of ribosome entry rate using ribosomal compartment 
#         t: time point (float)
#         kL: mitoribosome entry rate
#         kD: mitochondrial RNA Degradation rate
#     """
#     eps = 1e-16
#     if np.absolute(kL) <= eps: #limit kL = 0 case
#         Lambda = 1 - np.exp(-kD * t) * (1 + kD * t)
#     else:
#         Lambda = 1 - (kL + kD) / kL * np.exp(-kD * t) + kD / kL * np.exp(-(kD + kL) * t)
#     if Lambda < 0:#numerical under/overflow
#         Lambda = 0
#     return Lambda
