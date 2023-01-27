# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

from _my_functions import filename, pi_int

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

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

def g2_dressed_states(tau_in, Omega_in, gamma_in, w0_in):
    """
    Analytic expression for the normalised second-order correlation function
    of the dressed state transitions.
    """
    from numpy import exp, ones
    if abs(w0_in) == abs(Omega_in):
        g2_f = 1.0 - exp(-0.5 * gamma_in * tau_in)
    elif abs(w0_in) == 0.0:
        g2_f = ones(len(tau_in))
    return g2_f

#-----------------------------------------------------------------------------#
#                               FILENAME THINGS                               #
#-----------------------------------------------------------------------------#
# Read parameters
gamma, Omega, w0, kappa, dw, epsilon, N, phase = \
    np.genfromtxt(filename("parameters"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
Omega_str = pi_int(Omega)
w0_str = pi_int(w0)

# Pull data from file
# tau
tau = np.genfromtxt(filename("g2_corr"), usecols=0)
# g2
corr_filter = np.genfromtxt(filename("g2_corr"), usecols=1)

# analytic
corr_atom = g2_analytic(tau, Omega, gamma)
# Dressed states g2
corr_dressed = g2_dressed_states(tau, Omega, gamma, w0)

# print("Initial g2 value = {}".format(corr_filter[0]))

#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
plt.figure(figsize=[10, 6])

plt.plot(tau, corr_filter, color='C0', label='Filtered Correlation')
    
# Plot unfilterd correlation
plt.plot(tau, corr_atom, color='k', ls='dotted', alpha=0.5, label='Full Atom Correlation')

# Plot dressed_state correlation
plt.plot(tau, corr_dressed, color='k', ls='dashed', alpha=0.5, label='Dressed-State Approximation')

plt.xlabel(r'$\gamma \tau$', fontsize=15)
plt.ylabel(r'$g^{(2)}(\tau)$', fontsize=15)
plt.legend(loc='best', fontsize=15)
plt.title((r'$g^{{(2)}}_{{\mathrm{{Auto}}}}(\tau)$ with $\left( \Omega = {} \gamma,'
           r'N = {}, \delta\omega = {} \gamma, \kappa = {} \gamma, \omega_{{0}} = {}'
           r'\gamma \right)$').format(Omega_str, N, dw, kappa, w0_str),
          fontsize=15)

plt.tight_layout()
plt.show()
