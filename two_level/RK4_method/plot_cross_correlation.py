#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 15:15:50 2021

@author: jnga773
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
import sys
# Append path to add Jacobs_Function
sys.path.append("/home/jnga773/Documents/898-PhD/test/")
# sys.path.append("D:/Jacob/Google Drive/University/Physics/898 PhD/test")
from Jacobs_Functions import filename

plt.close('all')

def g2_analytic(tau_in, Omega_in, gamma_in):
    """
    Analytic expression for the normalised second-order correlation function
    from Howard 1
    """
    from numpy import exp, sqrt, cosh, sinh, real
    Omega_in = complex(Omega_in)
    gamma_in = complex(gamma_in)
    d_f = sqrt((0.25 * gamma_in) ** 2 - (Omega_in ** 2))
    g2_f = 1.0 - exp(-(3/4) * gamma_in * tau_in) * (cosh(d_f * tau_in) + \
                                                    ((3/4) * gamma_in / d_f) * sinh(d_f * tau_in))
    return real(g2_f)

#-------------------#
#     Read Data     #
#-------------------#

# Filenames
parent_dir = "./data_files/"
filename_param = filename("parameters", parent_dir)
filename_corr = filename("g2_corr", parent_dir)

# Pull parameters from file
gamma, Omega, w0a, w0b, kappaa, kappab, epsilon, N, dwa, dwb, phase, dt, t_max, tau1_max, tau2_max = \
    np.genfromtxt(fname=filename_param, delimiter="=", skip_header=1, usecols=(1,))
N = int(N)
phase = int(phase)
filters = (2*N) + 1
del(filename_param)

# Pull data from file
# tau
tau = np.genfromtxt(filename_corr, usecols=0)
# Auto-correlation of Cavity A
corr_autoA = np.genfromtxt(filename_corr, usecols=1)
# Auto-correlation of Cavity B
corr_autoB = np.genfromtxt(filename_corr, usecols=2)
# Cross correlation
corr_cross = np.genfromtxt(filename_corr, usecols=3)

#-------------------#
#     Plot Data     #
#-------------------#
# Plot titles
title_auto = r'Auto-Correlation $\left( \Omega = {} \gamma, N = {}, \epsilon = {} \right)$'
title_cross = (r'Cross-Correlation $\left( \omega_{{0}}^{{(A)}} = {}\gamma, '
              r'\kappa_{{a}} = {}\gamma, \delta\omega_{{A}} = {}\gamma, '
              r'\omega_{{0}}^{{(B)}} = {}\gamma, \kappa_{{B}} = {}\gamma, '
              r'\delta\omega_{{B}} = {}\gamma \right)$')
               
fig, ax = plt.subplots(2, 1, sharex=True, figsize=[12, 8])

# Auto-correlation
ax[0].plot(tau, corr_autoA, color='C0', ls='solid',
           label=r'$\langle A^{\dagger}(0) A^{\dagger} A(\tau) A(0) \rangle$')
ax[0].plot(tau, corr_autoB, color='C1', ls='dashed',
           label=r'$\langle B^{\dagger}(0) B^{\dagger} B(\tau) B(0) \rangle$')
ax[0].set_ylabel(r'$g^{(2)}_{\mathrm{auto}}(\tau)$', fontsize=15)
ax[0].set_title(title_auto.format(round(Omega, 3), N, epsilon), fontsize=15)
ax[0].legend(loc='best', fontsize=15)

# Cross-correlation
ax[1].plot(tau, corr_cross, color='C0', ls='solid',
           label=r'$\langle A^{\dagger}(0) B^{\dagger} B(\tau) A(0) \rangle$')
ax[1].plot(tau, 1 + np.exp(-0.5 * gamma * tau), color='k', ls='dashed', alpha=0.5,
           label='Dressed State Appoximation')
ax[1].set_ylabel(r'$g^{(2)}_{\mathrm{cross}}(\tau)$', fontsize=15)
ax[1].set_xlabel(r'$\gamma \tau$', fontsize=15)
ax[1].set_title(title_cross.format(round(w0a, 3), kappaa, dwa, round(w0b, 3), kappab, dwb), fontsize=15)
ax[1].legend(loc='best', fontsize=15)

fig.tight_layout()
fig.show()