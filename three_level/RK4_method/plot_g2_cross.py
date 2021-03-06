# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

from _my_functions import filename
from _dressed_state_correlations import three_level_eig

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

plt.close('all')

#-----------------------------------------------------------------------------#
#                               FILENAME THINGS                               #
#-----------------------------------------------------------------------------#
# Read parameters
gamma, Omega, alpha, delta, xi, w0a, kappaa, dwa, w0b, kappab, dwb, epsilon, N, phase = \
    np.genfromtxt(filename("cross_parameters"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
Omega = round(Omega, 2)
w0a = round(w0a, 2)
w0b = round(w0b, 2)

# Pull data from file
# tau
tau = np.genfromtxt(filename("cross_g2_corr"), usecols=0)
# g2
corr_auto = np.genfromtxt(filename("g2_corr"), usecols=1)
corr_cross = np.genfromtxt(filename("cross_g2_corr"), usecols=1)

#-----------------------------------------------------------------------------#
#                        PLOT AUTO AND CROSS g^{(2)}                          #
#-----------------------------------------------------------------------------#
fig, ax = plt.subplots(1, 2, figsize=[15, 6])

# Auto-correlation
ax[0].plot(tau, corr_auto, color='C0',
           label='Filtered Correlation $g^{{(2)}}(0) = {}$'.format(round(corr_auto[0], 2)))
# ax[0].plot(tau, dressed_auto, color='k', ls='dashed', alpha=0.5, label='Dressed State Approximation')
# ax[0].plot(tau, me_corr_auto, color='C1', ls='dashed', label='Master Equation')

ax[0].set_ylabel(r'$g^{(2)}_{\mathrm{Auto}}(\tau)$', fontsize=15)
ax[0].set_xlabel(r'$\gamma \tau$', fontsize=15)
ax[0].set_title((r'Auto-correlation with $N = {}, \delta\omega_{{a}} = {}'
                 r'\Gamma, \kappa_{{a}} = {} \Gamma, \omega_{{0}}^{{(a)}} = {}'
                 r'\Gamma$').format(N, dwa, kappaa, w0a), fontsize=15)
ax[0].legend(loc='best', fontsize=15)

# Cross-correlation
ax[1].plot(tau, corr_cross, color='C0',
           label='Filtered Correlation $g^{{(2)}}(0) = {}$'.format(round(corr_cross[0], 2)))
# ax[1].plot(tau, dressed_cross, color='k', ls='dashed', alpha=0.5, label='Dressed State Approximation')
# ax[1].plot(tau, me_corr_cross, color='C1', ls='dashed', label='Master Equation')

ax[1].set_ylabel(r'$g^{(2)}_{\mathrm{Cross}}(\tau)$', fontsize=15)
ax[1].set_xlabel(r'$\gamma \tau$', fontsize=15)
ax[1].set_title((r'Auto-correlation with $N = {}, \delta\omega_{{b}} = {}'
                 r'\Gamma, \kappa_{{b}} = {} \Gamma, \omega_{{0}}^{{(b)}} = {}'
                 r'\Gamma$').format(N, dwb, kappab, w0b), fontsize=15)
ax[1].legend(loc='best', fontsize=15)

fig.suptitle((r'Auto- and Cross-Correlations with $\Omega = {} \Gamma,'
              r'\alpha = {} \Gamma, \delta = {} \Gamma, \xi = {}$').format(Omega, alpha, delta, xi), fontsize=15)
fig.tight_layout()
fig.show()
