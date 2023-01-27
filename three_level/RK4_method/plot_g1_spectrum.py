# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

from Jacobs_Functions import filename, spectrum, norm_spectra
from analytic_functions._three_level_moments import g1_calc
from analytic_functions._dressed_state_correlations import three_level_eig, print_transition_frequencies, dressed_g1_calc

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

plt.close('all')

#------------------------------------------------------------------------------#
#                                FILENAME THINGS                               #
#------------------------------------------------------------------------------#
# Read parameters
Gamma, Omega, alpha, delta, xi, N, halfwidth, kappa, dw, w0 = \
    np.genfromtxt(filename("g1_parameters"), delimiter="=", usecols=1)
N = int(N)
# Omega = round(Omega, 2)
# w0 = round(w0, 2)
# xi = round(xi, 2)

# Print transition frequencies
print_transition_frequencies(Omega, alpha, delta, xi)

# Read data
tau = np.genfromtxt(filename("g1_corr"), usecols=0)
g1 = np.genfromtxt(filename("g1_corr"), usecols=1) + 1j * \
     np.genfromtxt(filename("g1_corr"), usecols=2)

# Atomic/unfiltered correlation
g1_atom = g1_calc(tau, Gamma, Omega, alpha, delta, xi)

# Dressed state correlation
g1_dressed = dressed_g1_calc(tau, w0, Gamma, Omega, alpha, delta, xi)

# Calculate Fourier transform
spec, wlist = spectrum(tau, g1, norm='peak')
spec_atom, wlist = spectrum(tau, g1_atom, norm='peak')
spec_dressed, wlist = spectrum(tau, g1_dressed, norm='peak')

# Renormalise filtered spectra
# spec = norm_spectra(spec_atom, spec, wlist)
# Renomalise dressed state spectra
# spec_dressed = norm_spectra(spec_atom, spec_dressed, wlist)

#-----------------------------------------------------------------------------#
#                               PLOT G^{(1)}                                  #
#-----------------------------------------------------------------------------#
fig, ax = plt.subplots(2, 1, sharex=True, figsize=[8, 6])

# Real part
ax[0].plot(tau, g1.real, color='C0', ls='solid', lw=2.0, label='Filtered')
ax[0].plot(tau, g1_dressed.real, color='C1', ls='dashed', lw=2.0, label='Dressed State')
# ax[0].plot(tau, g1_atom.real, color='k', ls='dashed', lw=1.0, label='Unfiltered')

ax[0].set_ylabel(r'Real Part', fontsize=15)
ax[0].set_xlim(-0.2, 5.2)
ax[0].legend(loc='best', fontsize=15)
ax[0].set_title((r'$G^{{(1)}}(\tau)$ with $\left( \Omega = {} \gamma,'
                 r'\alpha = {} \gamma, \delta = {} \gamma, \xi = {} \right)$'
                 ).format(Omega, alpha, delta, xi), fontsize=15)

ax[1].plot(tau, g1.imag, color='C0', ls='solid', lw=2.0)
ax[1].plot(tau, g1_dressed.imag, color='C1', ls='dashed', lw=2.0)
# ax[1].plot(tau, g1_atom.imag, color='k', ls='dashed', lw=1.0, alpha=0.5)

ax[1].set_ylabel(r'Imaginary Part', fontsize=15)
ax[1].set_xlabel(r'$\gamma t$', fontsize=15)
ax[1].set_title((r'$N = {}, \kappa = {} \gamma, \omega_{{0}} = {:.3f} \gamma,'
                 r'\delta\omega = {} \gamma$)'
                 ).format(N, kappa, w0, dw), fontsize=15)


fig.tight_layout()
fig.show()

#-----------------------------------------------------------------------------#
#                               PLOT SPECTRUM                                 #
#-----------------------------------------------------------------------------#
plt.figure(figsize=[8, 8])

plt.plot(wlist, spec, color='C0', ls='solid', lw=2.0,
         label='$N = {}, \delta\omega = {} \Gamma, \kappa = {} \Gamma$'.format(N, dw, kappa))
plt.plot(wlist, spec_dressed, color='C1', ls='dashed', lw=2.0, label='Dressed State')
plt.plot(wlist, spec_atom, color='k', ls='dashed', lw=1.0, label='Full Spectrum')

# Set lims for three-level or effective two-level
if 2 * delta == -alpha:
    plt.xlim(-1.5 * Omega, 1.5 * Omega)
else:
    plt.xlim(-abs(alpha), abs(alpha))
plt.xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=15)
plt.ylabel(r'Power Spectrum (a.u.)', fontsize=15)

plt.legend(loc='best', fontsize=15)
plt.title((r'$\Omega = {} \Gamma, \alpha = {} \Gamma, \delta = {} \Gamma, '
           r'\xi = {}, \omega_{{0}} = {:.3f} \Gamma$'
           ).format(Omega, alpha, delta, xi, w0), fontsize=15)

plt.tight_layout()
plt.show()

#-----------------------------------------------------------------------------#
#                          PLOT G^{(1)} AND SPECTRUM                          #
#-----------------------------------------------------------------------------#
# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 8])

# # First-order correlation
# ax[0].plot(tau, g1.real, color='C0', ls='solid', lw=2.0, label='Real')
# ax[0].plot(tau, g1.imag, color='C1', ls='dashed', lw=1.0, label='Imaginary')
# ax[0].set_xlim(-0.1, 10.1)

# ax[0].set_xlabel(r'$\gamma \tau$', fontsize=12)
# ax[0].set_ylabel(r'$G^{(1)}(\tau)$', fontsize=12)
# ax[0].set_title(r'$\Omega = {} \gamma, \alpha = {} \gamma, \delta = {} \gamma, \xi = {}$'.format(Omega, alpha, delta, xi), fontsize=12)
# ax[0].legend(loc='best', fontsize=12)

# # Spectra
# ax[1].plot(wlist, spec, color='C0', ls='solid', lw=2.0, label='Filtered')
# ax[1].plot(wlist, spec_me, color='k', ls='dashed', lw=1.0, alpha=0.5, label='Unfiltered')

# # Set lims for three-level or effective two-level
# if 2 * delta == -alpha:
#     ax[1].set_xlim(-1.5 * Omega, 1.5 * Omega)
# else:
#     ax[1].set_xlim(-abs(alpha), abs(alpha))
# ax[1].set_xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=12)
# ax[1].set_ylabel('Power Spectrum (a.u.)', fontsize=12)
# ax[1].set_title(r'$N = {}, \delta\omega = {} \gamma, \kappa = {} \gamma, \omega_{{0}} = {} \gamma$'.format(N, dw, kappa, w0), fontsize=12)
# ax[1].legend(loc='best', fontsize=12)

# fig.tight_layout()
# fig.show()

#-----------------------------------------------------------------------------#
#                        SPECTRUM AND ZOOMED SPECTRUM                         #
#-----------------------------------------------------------------------------#
# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 8])

# # Full-picture spectra
# ax[0].plot(wlist, spec, color='C0', ls='solid', lw=2.0, label='Filtered')
# ax[0].plot(wlist, spec_me, color='k', ls='dashed', lw=1.0, alpha=0.5, label='Unfiltered')

# ax[0].set_xlim(-abs(alpha), abs(alpha))
# ax[0].set_xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=12)
# ax[0].set_ylabel('Power Spectrum (a.u.)', fontsize=12)
# ax[0].set_title(r'$N = {}, \delta\omega = {} \gamma, \kappa = {} \gamma, \omega_{{0}} = {} \gamma$'.format(N, dw, kappa, w0), fontsize=12)
# ax[0].legend(loc='best', fontsize=12)


# # Zoomed-in spectra
# ax[1].plot(wlist, spec, color='C0', ls='solid', lw=2.0)
# ax[1].plot(wlist, spec_me, color='k', ls='dashed', lw=1.0, alpha=0.5)

# # Set lims for three-level or effective two-level
# if 2 * delta == -alpha:
#     ax[1].set_xlim(-1.5 * Omega, 1.5 * Omega)
# else:
#     ax[1].set_xlim(-abs(alpha), abs(alpha))
# ax[1].set_ylim(-0.001, 0.1)
# ax[1].set_xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=12)
# ax[1].set_ylabel('Power Spectrum (a.u.)', fontsize=12)
# ax[1].set_title(r'The same but zoomed in a bit', fontsize=12)


# fig.tight_layout()
# fig.show()
