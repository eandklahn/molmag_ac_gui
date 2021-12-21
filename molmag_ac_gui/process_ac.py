#std packages
import os
from subprocess import Popen, PIPE

#third-party packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import scipy.constants as sc
from scipy.optimize import curve_fit

from lmfit import Parameters, minimize as lmfit_minimize

"""
HELPER FUNCTIONS FOR PROCESSING AC-DATA
"""

def fit_Xp_Xpp_genDebye(v, Xp_data, Xpp_data):
    
    def genDebye_objective(params, v, data):
        
        Xp_data = data[0]
        Xpp_data = data[1]
        
        Xp_model = Xp_dataset(params, v)
        Xpp_model = Xpp_dataset(params, v)
        
        Xp_resid = Xp_data - Xp_model
        Xpp_resid = Xpp_data - Xpp_model
    
        return np.concatenate([Xp_resid, Xpp_resid])
            
    data = [Xp_data, Xpp_data]
    
    tau_init = (v[np.argmax(Xpp_data)]*2*np.pi)**-1
    fit_params = Parameters()
    fit_params.add('Xs', value=Xp_data[-1], min=0, max=np.inf)
    fit_params.add('Xt', value=Xp_data[0], min=0, max=np.inf)
    fit_params.add('tau', value=tau_init, min=0, max=np.inf)
    fit_params.add('alpha', value=0.1, min=0, max=1)
    
    out = lmfit_minimize(genDebye_objective, fit_params, args=(v, data))
    
    Xs = out.params['Xs'].value
    Xt = out.params['Xt'].value
    tau = out.params['tau'].value
    alpha = out.params['alpha'].value
    resid = out.residual.sum()
    
    if 'Could not estimate error-bars.' in out.message:
        # When fit errors are not available from experiment, set to arbitratrily small
        # value. The actual standard deviations will then be based on alpha.
        tau_fit_err = 0.0001*tau
    else:
        tau_idx = out.var_names.index('tau')
        tau_fit_err = np.sqrt(out.covar[tau_idx, tau_idx])
    
    return tau, tau_fit_err, alpha, Xs, Xt, resid

def tau_err_RC(tau, tau_fit_err, alpha, n_sigma=1):
    """
    Calculates the error in tau from the alpha-values as
    defined by Reta and Chilton (RC) in https://doi.org/10.1039/C9CP04301B
    """
    
    dtau = []
    for idx in range(len(tau)):
        t = tau[idx]
        te = tau_fit_err[idx]
        a = alpha[idx]
        
        dtau1 = np.abs(np.log10(t**-1)-np.log10((t*np.exp((n_sigma*1.82*np.sqrt(a))/(1-a)))**-1))
        dtau2 = np.abs(np.log10(t**-1)-np.log10((t*np.exp((-n_sigma*1.82*np.sqrt(a))/(1-a)))**-1))
        dtau3 = np.abs(np.log10(t**-1)-np.log10((t+te)**-1))
        dtau4 = np.abs(np.log10(t**-1)-np.log10((t-te)**-1))
        
        dtau.append(max([dtau1, dtau2, dtau3, dtau4]))
    
    return np.array(dtau)

def diamag_correction(H, H0, Mp, Mpp, m_sample, M_sample, Xd_sample, constant_terms=[], paired_terms=[]):
    """
    Calculates a diamagnetic correction of the data in Mp and Mpp and calculates
    the corresponding values of Xp and Xpp
    
    Input
    H: amplitude of AC field (unit: Oe)
    H0: strength of applied DC field (unit: Oe)
    Mp: in-phase magnetization (unit: emu)
    Mpp: out-of-phase magnetization (unit: emu)
    m_sample: sample mass in (unit: mg)
    M_sample: sample molar mass (unit: g/mol)
    Xd_sample: sample diamagnetic susceptibility in emu/(Oe*mol) from DOI: 10.1021/ed085p532
    constant_terms: terms to be subtracted directly from magnetization (unit: emu/Oe)
    paired_terms: list of tuples (tup) to be subtracted pairwise from magnetization
                  The terms must pairwise have the unit emu/Oe when multiplied,
                  fx. unit of tup[0] is emu/(Oe*<amount>) and unit of tup[1] is <amount>
    
    Output
    Mp_molar: molar in-phase magnetization, corrected for diamagnetic contribution (unit: emu/mol)
    Mpp_molar: molar out-of-phase magnetization, corrected for diamagnetic contribution (unit: emu/mol)
    Xp_molar: molar in-phase susceptibility, corrected for diamagnetic contribution (unit: emu/(Oe*mol))
    Xpp_molar: molar out-of-phase susceptibility, corrected for diamagnetic contribution (unit: emu/(Oe*mol))
    """
    # Old
    #Xp = (Mp - self.Xd_capsule*H - self.Xd_film*film_mass*H)*molar_mass/(sample_mass*H) - Xd_sample*molar_mass
    #Xpp = Mpp/(sample_mass*H)*molar_mass
    
    # NEW (documentation in docs with eklahn@chem.au.dk)
    # ---------------------------
    # Recalculate sample mass into g
    m_sample *= 10**-3
    # Calculate the molar amount of the sample
    n_sample = m_sample/M_sample
    
    sum_of_constants = sum(constant_terms)
    sum_of_pairwise = sum([tup[0]*tup[1] for tup in paired_terms])
    
    Mp_molar = (Mp - (sum_of_constants + sum_of_pairwise)*H - Xd_sample*H*n_sample)/n_sample
    Mpp_molar = Mpp/n_sample
    
    Xp_molar = Mp_molar/H
    Xpp_molar = Mpp_molar/H

    return Mp_molar, Mpp_molar, Xp_molar, Xpp_molar

def Xp_(v, Xs, Xt, tau, alpha):
    """
    Calculates the function X' [chi prime] as specified in Molecular Nanomagnets eq. 3.27
    This is the extended Debye model
    
    Input:
    v: frequency of the AC field
    Xs: adiabatic limit of susceptibility
    Xt: isothermal limit of susceptibility
    tau: relaxation time of the system
    alpha: width of relaxation time distribution
    
    Output
    Xp: the function value at the specified frequency
    """

    w = 2*np.pi*v
    
    Xp = Xs + (Xt - Xs)*(1 + (w*tau)**(1-alpha)*np.sin(np.pi*alpha/2))/(1 + 2*(w*tau)**(1-alpha)*np.sin(np.pi*alpha/2) + (w*tau)**(2-2*alpha))
    
    return Xp

def Xpp_(v, Xs, Xt, tau, alpha):
    """
    Calculates the function X'' [chi double-prime] as specified in Molecular Nanomagnets eq. 3.27
    This is the extended Debye model
    
    Input:
    v: frequency of the AC field
    Xs: adiabatic limit of susceptibility
    Xt: isothermal limit of susceptibility
    tau: relaxation time of the system
    alpha: width of relaxation time distribution
    
    Output
    Xpp: the function value at the specified frequency
    """

    w = 2*np.pi*v
    
    Xpp = (Xt - Xs)*(w*tau)**(1-alpha)*np.cos(np.pi*alpha/2)/(1 + 2*(w*tau)**(1-alpha)*np.sin(np.pi*alpha/2) + (w*tau)**(2-2*alpha))

    return Xpp

def Xp_dataset(params, v):
    """
    Wrapper for Xp_ to use with lmfit.
    
    Input:
    params: Instance of lmfit.Parameters containing elements Xs, Xt, alpha and tau as defined under function "Xp_"
    
    Output:
    See output of function "Xp_"
    """
    
    tau = params['tau']
    alpha = params['alpha']
    Xt = params['Xt']
    Xs = params['Xs']
    return Xp_(v, Xs, Xt, tau, alpha)
    
def Xpp_dataset(params, v):
    """
    Wrapper for Xpp_ to use with lmfit.
    
    Input:
    params: Instance of lmfit.Parameters containing elements Xs, Xt, alpha and tau as defined under function "Xpp_"
    
    Output:
    See output of function "Xpp_"
    """
    
    tau = params['tau']
    alpha = params['alpha']
    Xt = params['Xt']
    Xs = params['Xs']
    return Xpp_(v, Xs, Xt, tau, alpha)
    
def getParameterGuesses(T, tau):
    """
    DOCUMENTED
    Calculates guesses for optimal fitting parameters to begin the fit
    
    Input
    T: temperature array
    tau: relaxation time array
    
    Output
    guess: dictionary of guessed parameters for the relaxation functions
    """
    
    kB = sc.Boltzmann
    
    # Obtaining points for Orbach parameter guessing
    T1, tau1 = T[-1], tau[-1]
    T2, tau2 = T[-2], tau[-2]
    
    # Calculating guesses for Orbach parameters
    Ueff_guess = kB*np.log(tau2/tau1)/(1/T2-1/T1)
    
    t0_guess = tau1*np.exp(-Ueff_guess/(kB*T1))
    
    # Obtaining points for Raman parameter guessing
    l = len(T)
    if l/2 % 1 == 0:
        index1, index2 = int(l/2-1), int(l/2)
    else:
        index1, index2 = int(l/2-0.5), int(l/2+0.5)
    
    T1, tau1 = T[index1], tau[index1]
    T2, tau2 = T[index2], tau[index2]    
    
    # Calculating guesses for Raman parameters
    n_guess = (np.log(tau1) - np.log(tau2))/(np.log(T2) - np.log(T1))
    
    Cr_guess = (tau1**(-1))*(T1**(-n_guess))
    
    # Extracting guess for QT parameter
    tQT_guess = tau[0]
    
    guess = {'Ueff': Ueff_guess,
            't0': t0_guess,
            'n': n_guess,
            'Cr': Cr_guess,
            'tQT': tQT_guess}
    
    return guess

def _QT(T, tQT):
    """
    Basic function for calculating relaxation time due to
    quantum tunneling
    
    Input
    T: temperature for the calculation
    tQT: characteristic time for quantum tunneling
    
    Output
    tau: relaxation time due to quantum tunneling
    """
    
    tau = tQT
    
    return tau

def _R(T, Cr, n):
    """
    Basic function for calculating relaxation time due to
    the Raman mechanism. For canonical definition, see fx.
    DOI: 10.1039/c9cc02421b
    
    Input
    T: temperature for the calculation
    Cr: Raman pre-factor
    n: Raman exponent
    
    Output
    tau: relaxation time due to the Raman mechanism
    """
    
    tau = (Cr**(-1))*(T**(-n))

    return tau
    
def _O(T, t0, Ueff):
    """
    Basic function for calculating relaxation time due to
    the Orbach relaxation mechanism
    
    Input
    T: temperature for the calculation
    t0: characteristic time for quantum tunneling
    Ueff: effective energy barrier to thermal relaxation
    
    Output
    tau: relaxation time due to the Orbach mechanism
    """
    
    kB = sc.Boltzmann
    
    tau = t0*np.exp(Ueff/(kB*T))
    
    return tau
    
def _QT_log(T, tQT):
    """
    Wrapper function to function _QT that computes the logarithmic
    relaxation time due to quantum tunneling.
    
    See help(_QT) for more
    """
    
    return np.log(_QT(T, tQT))

def _R_log(T, Cr, n):
    """
    Wrapper function to function _R that computes the logarithmic
    relaxation time due to the Raman mechanism.
    
    See help(_R) for more
    """    
    
    return np.log(_R(T, Cr, n))
    
def _O_log(T, t0, Ueff):
    """
    Wrapper function to function _O that computes the logarithmic
    relaxation time due to the Orbach mechanism.
    
    See help(_O) for more
    """
    
    return np.log(_O(T, t0, Ueff))
    
def _QTR(T, tQT, Cr, n):
    """
    Wrapper function that computes the combined effect of a quantum
    tunneling mechanism and the Raman mechanism
    
    See help(_QT) and help(_R) for more
    """
    
    w = 1/_QT(T, tQT) + 1/_R(T, Cr, n)
    
    tau = 1/w
    
    return np.log(tau)

def _QTO(T, tQT, t0, Ueff):
    """
    Wrapper function that computes the combined effect of a quantum
    tunneling mechanism and the Orbach mechanism
    
    See help(_QT) and help(_O) for more
    """
    
    w = 1/_QT(T, tQT) + 1/_O(T, t0, Ueff)
    
    tau = 1/w
    
    return np.log(tau)
    
def _RO(T, Cr, n, t0, Ueff):
    """
    Wrapper function that computes the combined effect of a Raman
    mechanism and the Orbach mechanism
    
    See help(_R) and help(_O) for more
    """
    
    w = 1/_R(T, Cr, n) + 1/_O(T, t0, Ueff)
    
    tau = 1/w
    
    return np.log(tau)

def _QTRO(T, tQT, Cr, n, t0, Ueff):
    """
    Wrapper function that computes the combined effect of a quantum
    tunneling mechanism, the Raman mechanism and the Orbach mechanism
    
    See help(_QT), help(_R) and help(_O) for more
    """
    
    w = 1/_QT(T, tQT) + 1/_R(T, Cr, n) + 1/_O(T, t0, Ueff)
    
    tau = 1/w
    
    return np.log(tau)

def getStartParams(guess, fitType='QTRO'):
    
    p0 = 0
    if fitType=='QT':
        p0 = [guess['tQT']]
    elif fitType=='R':
        p0 = [guess['Cr'], guess['n']]
    elif fitType=='O':
        p0 = [guess['t0'], guess['Ueff']]
    elif fitType=='QTR':
        p0 = getStartParams(guess, fitType='QT') + getStartParams(guess, fitType='R')
    elif fitType=='QTO':
        p0 = getStartParams(guess, fitType='QT') + getStartParams(guess, fitType='O')
    elif fitType=='RO':
        p0 = getStartParams(guess, fitType='R') + getStartParams(guess, fitType='O')
    elif fitType=='QTRO':
        p0 = [guess['tQT'], guess['Cr'], guess['n'], guess['t0'], guess['Ueff']]
    else:
        print('fitType parameter did not correspond to any correct one')
        
    return p0

def getFittingFunction(fitType='QTRO'):
    
    f = 0
    if fitType=='QT':
        f = _QT_log
    elif fitType=='R':
        f = _R_log
    elif fitType=='O':
        f = _O_log
    elif fitType=='QTR':
        f = _QTR
    elif fitType=='QTO':
        f = _QTO
    elif fitType=='RO':
        f = _RO
    elif fitType=='QTRO':
        f = _QTRO
    else:
        print('fitType parameter did not correspond to any correct one')
        
    return f
    
def readPopt(popt, pcov, fitType='QTRO'):

    p_fit = 0
    if fitType=='QT':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT']}
    elif fitType=='R':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['Cr', 'n']}
    elif fitType=='O':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['t0', 'Ueff']}
    elif fitType=='QTR':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT', 'Cr', 'n']}
    elif fitType=='QTO':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT', 't0', 'Ueff']}
    elif fitType=='RO':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['Cr', 'n', 't0', 'Ueff']}
    elif fitType=='QTRO':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT', 'Cr', 'n', 't0', 'Ueff']}
    
    return p_fit

def readPFITinOrder(p_fit, plotType='O'):

    p = []
    if plotType=='QT':
        tQT = p_fit['params'][p_fit['quantities'].index('tQT')]
        p = [tQT]
    elif plotType=='R':
        Cr = p_fit['params'][p_fit['quantities'].index('Cr')]
        n = p_fit['params'][p_fit['quantities'].index('n')]
        p = [Cr, n]
    elif plotType=='O':
        t0 = p_fit['params'][p_fit['quantities'].index('t0')]
        Ueff = p_fit['params'][p_fit['quantities'].index('Ueff')]
        p = [t0, Ueff]
    elif plotType=='QTR':
        p = readPFITinOrder(p_fit, plotType='QT') + readPFITinOrder(p_fit, plotType='R')
    elif plotType=='QTO':
        p = readPFITinOrder(p_fit, plotType='QT') + readPFITinOrder(p_fit, plotType='O')
    elif plotType=='RO':
        p = readPFITinOrder(p_fit, plotType='R') + readPFITinOrder(p_fit, plotType='O')
    elif plotType=='QTRO':
        p = readPFITinOrder(p_fit, plotType='QT') + readPFITinOrder(p_fit, plotType='R') + readPFITinOrder(p_fit, plotType='O')
    
    return p
    
def addPartialModel(fig, Tmin, Tmax, p_fit, plotType='O', *args, **kwargs):
    
    ax = fig.get_axes()[0]
    
    f = getFittingFunction(fitType=plotType)
    p = readPFITinOrder(p_fit, plotType=plotType)

    #Quickfix for zero division when calculating: line = ax.plot(1/T_space ...) 
    if Tmin == 0: 
        Tmin = 0.00001 
    
    T_space = np.linspace(Tmin, Tmax, 101)
    line, = ax.plot(1/T_space, np.ones(T_space.shape)*f(T_space, *p), *args, **kwargs)
    return line
    