#std packages
import os
from subprocess import Popen, PIPE

#third-party packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import scipy.constants as sc
import warnings


from lmfit import Parameters, minimize as lmfit_minimize

kB = sc.Boltzmann

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
    warnings.filterwarnings("ignore", category=RuntimeWarning, module=__name__)

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
    Htot = H + H0
    
    m_sample *= 10**-3
    # Calculate the molar amount of the sample
    n_sample = m_sample/M_sample
    
    sum_of_constants = sum(constant_terms)
    sum_of_pairwise = sum([tup[0]*tup[1] for tup in paired_terms])
    
    Mp_molar = (Mp - (sum_of_constants + sum_of_pairwise)*Htot - Xd_sample*Htot*n_sample)/n_sample
    Mpp_molar = Mpp/n_sample
    
    Xp_molar = Mp_molar/Htot
    Xpp_molar = Mpp_molar/Htot

    return Mp_molar, Mpp_molar, Xp_molar, Xpp_molar


def diamag_correction_dc(H, H0, M, m_sample, M_sample):
    """
    This should be fixed. 
    Calculates a diamagnetic correction of the data in Mp and Mpp and calculates
    the corresponding values of Xp and Xpp
    
    Input
    H0: strength of applied DC field (unit: Oe)
    M: Moment (unit: emu)
    m_sample: sample mass in (unit: mg)
    M_sample: sample molar mass (unit: g/mol)
    Xd_sample: sample diamagnetic susceptibility in emu/(Oe*mol) from DOI: 10.1021/ed085p532

    Output
    M_molar: molar moment (unit: emu/mol)
    X_molar: molar in-phase susceptibility (unit: emu/(Oe*mol)). I am unsure if this should included 
    """
    Htot = H + H0
    Xd_sample = 0 #THIS SHOULD BE CHANGED / RECONSIDERED 
    m_sample *= 10**-3
    n_sample = m_sample/M_sample
    M_molar = M/n_sample
    X_molar = M_molar/Htot - Xd_sample

    return M_molar, X_molar


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
    
def getParameterGuesses(T, tau, fitwith):
    """
    DOCUMENTED
    Calculates guesses for optimal fitting parameters to begin the fit
    
    Input
    T: temperature array
    tau: relaxation time array
    
    Output
    guess: dictionary of guessed parameters for the relaxation functions
    """
    scale = 10**20
    kB = sc.Boltzmann
    
    # Obtaining points for Orbach parameter guessing
    T1, tau1 = T[-1], tau[-1]
    T2, tau2 = T[-2], tau[-2]
    
    # Calculating guesses for Orbach parameters
    Ueff_by_kB_guess = np.log(tau2/tau1)/(1/T2-1/T1)
    
    t0_scaled_guess = tau1*10**7*np.exp(-Ueff_by_kB_guess/T1)
    
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
    
    if fitwith == 'all':
        fitwith = 'QTRO'

    params = Parameters()
    params.add(name='tQT', value=tQT_guess, vary='QT' in fitwith, min=0)
    params.add(name='Cr', value=Cr_guess, vary='R' in fitwith, min=0)
    params.add(name='n', value=n_guess, vary='R' in fitwith, min=0, max=20)
    params.add(name='t0_scaled', value=t0_scaled_guess, vary='O' in fitwith,  min = 10**-3, max = 10**6)
    params.add(name='Ueff_by_kB', value=Ueff_by_kB_guess, vary='O' in fitwith, min = 10**-3, max = 1000)
    params.add(name='useQT', value=int('QT' in fitwith), vary=False, min = 0, max = 1)
    params.add(name='useR', value=int('R' in fitwith), vary=False, min = 0, max = 1)
    params.add(name='useO', value=int('O' in fitwith), vary=False, min = 0, max = 1)

    return params

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
    
def _O(T, t0_scaled, Ueff_by_kB):
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

    tau = t0_scaled*10**-7*np.exp(Ueff_by_kB/T)
    
    return tau
    

def getStartParams(guess, fitType='QTRO'):
    
    p0 = 0
    if fitType=='QT':
        p0 = [guess['tQT']]
    elif fitType=='R':
        p0 = [guess['Cr'], guess['n']]
    elif fitType=='O':
        p0 = [guess['t0_scaled'], guess['Ueff_by_kB']]
    elif fitType=='QTR':
        p0 = getStartParams(guess, fitType='QT') + getStartParams(guess, fitType='R')
    elif fitType=='QTO':
        p0 = getStartParams(guess, fitType='QT') + getStartParams(guess, fitType='O')
    elif fitType=='RO':
        p0 = getStartParams(guess, fitType='R') + getStartParams(guess, fitType='O')
    elif fitType=='QTRO':
        p0 = [guess['tQT'], guess['Cr'], guess['n'], guess['t0_scaled'], guess['Ueff_by_kB']]
    else:
        print('fitType parameter did not correspond to any correct one')
        
    return p0


    
def readPopt(popt, pcov, fitType='QTRO'):

    p_fit = 0
    if fitType=='QT':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT']}
    elif fitType=='R':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['Cr', 'n']}
    elif fitType=='O':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['t0_scaled', 'Ueff_by_kB']}
    elif fitType=='QTR':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT', 'Cr', 'n']}
    elif fitType=='QTO':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT', 't0_scaled', 'Ueff_by_kB']}
    elif fitType=='RO':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['Cr', 'n', 't0_scaled', 'Ueff_by_kB']}
    elif fitType=='QTRO':
        p_fit = {'params': popt, 'sigmas': np.sqrt(np.diag(pcov)), 'quantities': ['tQT', 'Cr', 'n', 't0_scaled', 'Ueff_by_kB']}
    
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
        t0_scaled = p_fit['params'][p_fit['quantities'].index('t0_scaled')]
        Ueff_by_kB = p_fit['params'][p_fit['quantities'].index('Ueff_by_kB')]
        p = [t0_scaled, Ueff_by_kB]
    elif plotType=='QTR':
        p = readPFITinOrder(p_fit, plotType='QT') + readPFITinOrder(p_fit, plotType='R')
    elif plotType=='QTO':
        p = readPFITinOrder(p_fit, plotType='QT') + readPFITinOrder(p_fit, plotType='O')
    elif plotType=='RO':
        p = readPFITinOrder(p_fit, plotType='R') + readPFITinOrder(p_fit, plotType='O')
    elif plotType=='QTRO':
        p = readPFITinOrder(p_fit, plotType='QT') + readPFITinOrder(p_fit, plotType='R') + readPFITinOrder(p_fit, plotType='O')
    
    return p
    


def add_partial_model(fig, Tmin, Tmax, params, N_points=1000, *args, **kwargs):
    """
    Reimplementation of addPartialModel that supports the use of
    lmfit.Parameters for setting the function to plot.
    """

    ax = fig.get_axes()[0]
    T = np.linspace(Tmin, Tmax, N_points)
    #Calculates tau based on the parameters
    tau = relaxation_dataset(params, T)

    line, = ax.plot(1/T, np.log(tau), *args, **kwargs)

    return line

def general_relaxation(T, tQT, Cr, n, t0_scaled, Ueff_by_kB, useQT, useR, useO):

    tauQT = _QT(T, tQT)
    tauO = _O(T, t0_scaled, Ueff_by_kB)
    tauR = _R(T, Cr, n)

    rate = useQT*(1/tauQT)+useO*(1/tauO)+useR*(1/tauR) # We can add rates, but not relaxation times.
    tau = 1/rate 

    return tau

def relaxation_dataset(params, T):
    tQT = params['tQT']
    Cr = params['Cr']
    n = params['n']
    t0_scaled = params['t0_scaled']
    Ueff_by_kB = params['Ueff_by_kB']
    useQT = params['useQT']
    useR = params['useR']
    useO = params['useO']

    return general_relaxation(T, tQT, Cr, n, t0_scaled, Ueff_by_kB, useQT, useR, useO)

def fit_relaxation(T, tau_data, params):

    def objective(params, T, tau_data):
        residual = np.log(tau_data) - np.log(relaxation_dataset(params, T))
        return residual
    

    out = lmfit_minimize(objective, params, args=(T, tau_data), method = 'least_squares')

    return out

def default_parameters(fitwith='QTRO'):

    params = Parameters()
    params.add(name='tQT', value=1e-4, vary='QT' in fitwith, min=0)
    params.add(name='Cr', value=1e-3, vary='R' in fitwith, min=0)
    params.add(name='n', value=4, vary='R' in fitwith, min=0, max=20)
    params.add(name='t0_scaled', value=1, vary='O' in fitwith, min=0.001, max = 10000)
    params.add(name='Ueff_by_kB', value=50, vary='O' in fitwith, min = 10**-3, max = 1000)
    params.add(name='useQT', value=int('QT' in fitwith), vary=False, min = 0, max = 1)
    params.add(name='useR', value=int('R' in fitwith), vary=False, min = 0, max = 1)
    params.add(name='useO', value=int('O' in fitwith), vary=False, min = 0, max = 1)

    return params