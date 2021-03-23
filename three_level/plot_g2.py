# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
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

#-----------------------------------------------------------------------------#
#                               FILENAME THINGS                               #
#-----------------------------------------------------------------------------#
# Read parameters
gamma, Omega, alpha, delta, xi, w0, kappa, dw, epsilon, N, phase = \
    np.genfromtxt(filename("parameters"), delimiter="=", skip_header=1, usecols=1)
xi = round(xi, 2)
N = int(N)
phase = int(phase)
Omega = round(Omega, 2)
w0 = round(w0, 2)

# Read the data
tau = np.genfromtxt(filename("g2_corr"), usecols=0)
g2 = np.genfromtxt(filename("g2_corr"), usecols=1)

# Master equation to compare
# g2_me = np.genfromtxt(filename("me_g2_corr"), usecols=1)

# Unfiltered correlation
g2_atom = np.genfromtxt(filename("atom_g2_corr"), usecols=1)

#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
plt.figure(figsize=[10, 6])

plt.plot(tau, g2, color='C0', label='Filtered Correlation')
# plt.plot(tau, g2_me, color='C1', ls='--', label='Master Equation')

# Plot unfilterd correlation
plt.plot(tau, g2_atom, color='k', alpha=0.5, ls='dashed', label='Unfiltered Correlation')

# Plot dressed_state correlation
# plt.plot(tau, corr_dressed, ls='--', color='k', alpha=0.5, label='Dressed-State Approximation')

plt.xlabel(r'$\gamma \tau$', fontsize=15)
plt.ylabel(r'$g^{(2)}(\tau)$', fontsize=15)
plt.legend(loc='best', fontsize=15)
plt.title(r'$\left( \Omega = {}, \alpha = {}, \delta = {}, \xi = {} \right), \left( \kappa = {}, \omega_{{0}} = {}, N = {}, \delta\omega = {} \right)$'.format(Omega, alpha, delta, xi, kappa, w0 ,N, dw), fontsize=15)
plt.tight_layout()
plt.show()
