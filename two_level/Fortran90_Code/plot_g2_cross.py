# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

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

def g2_dressed_states(tau_in, Omega_in, gamma_in, w0_in, auto_or_cross='auto'):
    """
    Analytic expression for the normalised second-order correlation function
    of the dressed state transitions.
    """
    from numpy import exp, ones
    # If auto-correlation, calculate straight or anti-bunched
    if auto_or_cross == 'auto':
        if abs(w0_in) == abs(Omega_in):
            g2_f = 1.0 - exp(-0.5 * gamma_in * tau_in)
        elif abs(w0_in) == 0.0:
            g2_f = ones(len(tau_in))
    elif auto_or_cross == 'cross':
        g2_f = 1 + exp(-0.5 * gamma_in * tau_in)
    return g2_f

#-----------------------------------------------------------------------------#
#                               FILENAME THINGS                               #
#-----------------------------------------------------------------------------#
# Read parameters
gamma, Omega, N, halfwidth, kappa, dw, w0a, w0b = \
    np.genfromtxt(filename("g2_cross_parameters"), delimiter="=", usecols=1)
N = int(N)
phase = 1
# Omega = round(Omega, 2)
# w0a = round(w0a, 2)
# w0b = round(w0b, 2)

if Omega == 5 * np.pi:
    Omega_str = r'5 \pi'
if w0a == 5 * np.pi:
    w0a_str = r'5 \pi'
if w0b == -5 * np.pi:
    w0b_str = r'-5 \pi'

# Pull data from file
# tau
tau = np.genfromtxt(filename("g2_cross_corr"), usecols=0)
# g2
# corr_auto = np.genfromtxt(filename("g2_auto_corr"), usecols=1)
corr_cross = np.genfromtxt(filename("g2_cross_corr"), usecols=1)

# analytic
corr_atom = g2_analytic(tau, Omega, gamma)
# Dressed states g2
# dressed_auto = g2_dressed_states(tau, Omega, gamma, w0a, 'auto')
dressed_cross = g2_dressed_states(tau, Omega, gamma, w0a, 'cross')

# print("Initial g2 value = {}".format(corr_filter[0]))

#-----------------------------------------------------------------------------#
#                              PLOT CROSS g^{(2)}                             #
#-----------------------------------------------------------------------------#
# Plot
plt.figure(figsize=[8, 6])

if N > 0:
    plt.plot(tau, corr_cross, color='C0', ls='solid',
             label=(r'Multi-Mode $\left( N = {}, \delta\omega = {} \gamma,'
                    r'\kappa = {} \gamma \right)$').format(N, dw, kappa))
    plt.title((r'$g^{{(2)}}_{{\mathrm{{Cross}}}}(\tau)$ with $\left( \Omega = {}'
               r'\gamma, \omega_{{0}}^{{(a)}} = {} \gamma, \omega_{{0}}^{{(b)}}'
               r'= {} \gamma \right)$'
               ).format(Omega_str, w0a_str, w0b_str),
              fontsize=15)
if N == 0:
    plt.plot(tau, corr_cross, color='C0', ls='solid',
             label=(r'Single-Mode $\left( \kappa = {} \gamma \right)$').format(kappa))
    plt.title((r'$g^{{(2)}}_{{\mathrm{{Cross}}}}(\tau)$ with'
               r'$\left( \Omega = {} \gamma, \omega_{{0}}^{{(a)}} = {} \gamma,'
               r'\omega_{{0}}^{{(b)}} = {} \gamma \right) $'
               ).format(Omega_str, w0a_str, w0b_str),
              fontsize=15)
    
plt.plot(tau, dressed_cross, color='k', ls='dashed', alpha=0.5, label='Dressed State Approximation')

plt.xlabel(r'$\gamma \tau$', fontsize=15)
plt.ylabel(r'$g^{(2)}_{\mathrm{cross}}(\tau)$', fontsize=15)
plt.legend(loc='upper right', fontsize=15)

plt.tight_layout()
plt.show()

#-----------------------------------------------------------------------------#
#                         PLOT AUTO AND CROSS g^{(2)}                         #
#-----------------------------------------------------------------------------#
# fig, ax = plt.subplots(1, 2, figsize=[15, 6])

# # Auto-correlation
# ax[0].plot(tau, corr_auto, color='C0', label='Filtered Correlation')
# ax[0].plot(tau, dressed_auto, color='k', ls='dashed', alpha=0.5, label='Dressed State Approximation')

# ax[0].set_ylabel(r'$g^{(2)}_{\mathrm{Auto}}(\tau)$', fontsize=12)
# ax[0].set_xlabel(r'$\gamma \tau$', fontsize=15)
# ax[0].set_title(r'Auto-correlation with $\left( N = {}, \delta\omega_{{a}} = {}, \kappa_{{a}} = {} \gamma, \omega_{{0}}^{{(a)}} = {} \gamma \right)$'.format(N, dwa, kappaa, w0a), fontsize=15)
# ax[0].legend(loc='best', fontsize=15)

# # Crosscorrelation
# ax[1].plot(tau, corr_cross, color='C0', label='Filtered Correlation')
# ax[1].plot(tau, dressed_cross, color='k', ls='dashed', alpha=0.5, label='Dressed State Approximation')

# ax[1].set_ylabel(r'$g^{(2)}_{\mathrm{Cross}}(\tau)$', fontsize=12)
# ax[1].set_xlabel(r'$\gamma \tau$', fontsize=15)
# ax[1].set_title(r'Auto-correlation with $\left( N = {}, \delta\omega_{{b}} = {}, \kappa_{{b}} = {} \gamma, \omega_{{0}}^{{(b)}} = {} \gamma \right)$'.format(N, dwb, kappab, w0b), fontsize=15)
# ax[1].legend(loc='best', fontsize=15)

# fig.suptitle(r'Auto- and Cross-Correlations with $\Omega = {} \gamma$'.format(Omega), fontsize=15)
# fig.tight_layout()
# fig.show()

#-----------------------------------------------------------------------------#
#                          PLOT CROSS g^{(2)} COMPARE                         #
#-----------------------------------------------------------------------------#
# OG_data = np.genfromtxt(filename("cross_g2_corr_OG"), usecols=1)

# # Plot
# plt.figure(figsize=[6, 6])

# plt.plot(tau, OG_data, color='C0', ls='solid', label='OG Correct Data')
# plt.plot(tau, corr_cross, color='C1', ls='dashed', label='Newer Version')

# plt.xlabel(r'$\gamma \tau$', fontsize=12)
# plt.ylabel(r'$g^{(2)}_{\mathrm{cross}}(\tau)$', fontsize=12)
# plt.title(r'Comparing Program Versions', fontsize=12)

# plt.tight_layout()
# plt.show()

#-----------------------------------------------------------------------------#
#                        PLOT TWO-SIDED CROSS g^{(2)}                         #
#-----------------------------------------------------------------------------#
# # Flip and combine the cross-correlation data
# corr_two_sided = np.concatenate((np.flip(corr_cross[1:]), corr_cross))
# tau_two_sided = np.concatenate((-np.flip(tau[1:]), tau))

# # Plot
# plt.figure(figsize=[8, 6])
# plt.plot(tau_two_sided, corr_two_sided)

# plt.xlabel(r'$\gamma \tau$', fontsize=12)
# plt.ylabel(r'$g^{(2)}_{\mathrm{cross}}(\tau)$', fontsize=12)
# # plt.title(r'Comparing Program Versions', fontsize=12)

# plt.tight_layout()
# plt.show()