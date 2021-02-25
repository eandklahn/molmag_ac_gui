import numpy as np
import os
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scipy.constants as sc
from scipy.optimize import curve_fit

globalPointMarker = 'o'
globalMarkerSize = 5
globalTextSize = 20

def calcTcolor(T, Tmin, Tmax):
    """
    Calculates the color that the corresponding line should be plotted with based on the
    idea that the coldest temperature should be completely blue RGB=(0,0,255) and the warmest
    temperature should be completely red RGB=(255,0,0)
    
    Input
    T: temperature of the curve that is being plotted
    D: dictionary containing all temperatures

    Output
    RGB: a tuple containing (R,G,B)-values for the color to be used
    """
    
    T18 = Tmin + 1/8*(Tmax-Tmin)
    T38 = Tmin + 3/8*(Tmax-Tmin)
    T58 = Tmin + 5/8*(Tmax-Tmin)
    T78 = Tmin + 7/8*(Tmax-Tmin)
    
    R, B, G = 0, 0, 0
    if T >= T78:
        p = (T-T78)/(Tmax-T78)
        R = 1-0.5*p
    elif T >= T58:
        p = (T-T58)/(T78-T58)
        R = 1
        G = 1-p
    elif T >= T38:
        p = (T-T38)/(T58-T38)
        R = p
        G = 1
        B = 1-p
    elif T >= T18:
        p = (T-T18)/(T38-T18)
        G = p
        B = 1
    else:
        p = (T-Tmin)/(T18-Tmin)
        B = 0.5+0.5*p
    
    print('Here is the value of ', p)
    
    return (R,G,B)  
    
def Xp_(v, Xs, Xt, tau, alpha):
    """
    Calculates the function X' [chi prime] as specified in Molecular Nanomagnets eq. 3.27
    
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
    Ueff_guess = (np.log(tau2) - np.log(tau1))/(1/T2-1/T1)*kB
    
    t0_guess = tau1*np.exp(-Ueff_guess/(kB*T1))
    
    # Obtaining points for Raman parameter guessing
    l = len(T)
    index1, index2 = 0,0
    if l/2 % 1 == 0:
        index1, index2 = int(l/2), int(l/2+1)
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
    
def addPartialModel(fig, Tmin, Tmax, p_fit, plotType='O'):
    
    ax = fig.get_axes()[0]
    
    f = getFittingFunction(fitType=plotType)
    p = readPFITinOrder(p_fit, plotType=plotType)
    
    T_space = np.linspace(Tmin, Tmax, 101)
    line, = ax.plot(1/T_space, np.ones(T_space.shape)*f(T_space, *p), label=plotType)
    
    return line