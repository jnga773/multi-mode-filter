# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

import numpy as np
import matplotlib.pyplot as plt

from Jacobs_Functions import filename
from analytic_functions._dressed_state_correlations import print_transition_frequencies, dressed_g2_calc

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.close('all')

#-----------------------------------------------------------------------------#
#                               FILENAME THINGS                               #
#-----------------------------------------------------------------------------#
# Read parameters
Gamma, Omega, alpha, delta, xi, N, halfwidth, kappa, dw, w0 = \
    np.genfromtxt(filename("g2_auto_parameters"), delimiter="=", usecols=1)
xi = round(xi, 2)
N = int(N)

# Read the data
tau = np.genfromtxt(filename("g2_auto_corr"), usecols=0)
g2 = np.genfromtxt(filename("g2_auto_corr"), usecols=1)

# Print transition frequencies
print_transition_frequencies(Omega, alpha, delta, xi)

# Dressed-state correlations
g2_dressed = dressed_g2_calc(tau, w0, Gamma, Omega, alpha, delta, xi)

#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
plt.figure(figsize=[10, 6])

plt.plot(tau, g2, color='C0',
         label=(r'Filtered Correlation $\left( N = {}, \delta\omega = {}'
                r'\Gamma, \kappa = {} \Gamma, \omega_{{0}} = {} \Gamma \right)$'
                ).format(N, round(dw, 2), kappa, round(w0, 2)))

# Plot dressed state correlation
plt.plot(tau, g2_dressed, color='k', alpha=0.5, ls='dashed', label='Dressed State Correlation')

plt.xlabel(r'$\Gamma \tau$', fontsize=15)
plt.ylabel(r'$g^{(2)}(\tau)$', fontsize=15)

plt.title((r'$\left( \Omega = {} \Gamma, \alpha = {}, \delta = {}, \xi = {}'
           r'\right)$').format(round(Omega, 2), round(alpha,2), round(delta, 2), xi),
          fontsize=15)

plt.legend(loc='best', fontsize=15)

plt.tight_layout()
plt.show()
