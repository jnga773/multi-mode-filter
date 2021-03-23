#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 11:04:55 2021

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
from Jacobs_Functions import filename, spectrum

plt.close('all')

# Three-level Hamiltonian eigenvalues
def three_level_eigval(Omega_in, alpha_in, delta_in, xi_in):
    """
    Calculates the dressed/eigenvalues of the Hamiltonian in matrix form.
    """
    from numpy.linalg import eigvals
    from numpy import matrix, sort
    # Set up matrix
    H = matrix([[0, 0.5 * Omega_in, 0],
                [0.5 * Omega_in, -(0.5 * alpha_in) - delta_in, 0.5 * xi_in * Omega_in],
                [0, 0.5 * xi_in * Omega_in, -2 * delta_in]])
    # Calculate eigenvalues
    eig_out = eigvals(H)
    # Order them from big to small
    eig_out = sort(eig_out)
    return eig_out

# Normalisation function
def norm_spectra(Omega_in, alpha_in, delta_in, xi_in,
                 spec_atom_in, spec_in, w0_in):
    """
    Depending on where the filter is centred, normalise the filtered spectrum
    to the relevant peak.
    """
    from numpy import array
    # Calculate where the peaks are located in the UNfiltered spectrum
    from scipy.signal import find_peaks

    # Calculate eigenvalues of Three-Level Hamiltonian
    if 2 * delta_in == -alpha_in:
        peaks = find_peaks(spec_atom_in, height=0.1)[1]["peak_heights"]
        # Effectively a two-level atom so find the three-Mollow peaks
        peaks_calculated = [-Omega_in, 0.0, Omega_in]
    else:
        peaks = find_peaks(spec_atom_in, height=0.01)[1]["peak_heights"]
        # Usual three-level atom
        wm, w0, wp = three_level_eigval(Omega_in, alpha_in, delta_in, xi_in)

        # Calculate where the peaks should fall from eigenvalues
        peaks_calculated = array([-(wp - wm), -wp, wm, 0, -wm, wp, (wp - wm)])
        print("Calculated eigenvalue peaks are: ", peaks_calculated)

    # Cycle through peaks_calculated and compare with w0
    norm_peak_value = 1.0
    for place, item in enumerate(peaks_calculated):
        if round(item, 2) == round(w0_in, 2):
            # save the peak value
            norm_peak_value = peaks[place]
            # Leave loop
            break
    # Renormalise spectra
    spec_out = spec_in * norm_peak_value
    return spec_out

#------------------------------------------------------------------------------#
#                                FILENAME THINGS                               #
#------------------------------------------------------------------------------#
# Read parameters
gamma, Omega, alpha, delta, xi, w0, kappa, dw, epsilon, N, phase = \
    np.genfromtxt(filename("parameters"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
Omega = round(Omega, 2)
w0 = round(w0, 2)
xi = round(xi, 2)

# Read first-order correlation data
tau1 = np.genfromtxt(filename("g1_corr"), usecols=0)
g1 = np.genfromtxt(filename("g1_corr"), usecols=1) + 1j * \
     np.genfromtxt(filename("g1_corr"), usecols=2)
# Atomic/unfiltered first-order correlation
g1_atom = np.genfromtxt(filename("atom_g1_corr"), usecols=1) + 1j * \
          np.genfromtxt(filename("atom_g1_corr"), usecols=2)

# Calculate Fourier transform
spec, wlist = spectrum(tau1, g1, norm='peak')
spec_atom, wlist = spectrum(tau1, g1_atom, norm='peak')

# Renormalise filtered spectra
spec = norm_spectra(Omega, alpha, delta, xi, spec_atom, spec, w0)

# Read second-order correlation data
tau2 = np.genfromtxt(filename("g2_corr"), usecols=0)
g2 = np.genfromtxt(filename("g2_corr"), usecols=1)
# Unfiltered correlation
g2_atom = np.genfromtxt(filename("atom_g2_corr"), usecols=1)

# -----------------------------------------------------------------------------#
#                           PLOT g^{(1)} AND SPECTRUM                          #
# -----------------------------------------------------------------------------#
# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 8])

# # First-order correlation
# ax[0].plot(tau1, g1.real, color='C0', ls='solid', lw=2.0, label='Real')
# ax[0].plot(tau1, g1.imag, color='C1', ls='dashed', lw=1.0, label='Imaginary')
# ax[0].set_xlim(-0.1, 10.1)

# ax[0].set_xlabel(r'$\gamma \tau$', fontsize=12)
# ax[0].set_ylabel(r'$G^{(1)}(\tau)$', fontsize=12)
# ax[0].set_title(r'$\Omega = {} \gamma, \alpha = {} \gamma, \delta = {} \gamma, \xi = {}$'.format(Omega, alpha, delta, xi), fontsize=12)
# ax[0].legend(loc='best', fontsize=12)

# # Spectra
# ax[1].plot(wlist, spec, color='C0', ls='solid', lw=2.0, label='Filtered')
# ax[1].plot(wlist, spec_atom, color='k', ls='dashed', lw=1.0, alpha=0.5, label='Unfiltered')

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
#                               PLOT SPECTRUM                                 #
#-----------------------------------------------------------------------------#
# plt.figure(figsize=[8, 8])

# plt.plot(wlist, spec, color='C0', ls='solid', lw=2.0, label='Filtered Spectrum')
# plt.plot(wlist, spec_atom, color='k', ls='dashed', lw=1.0, alpha=0.5, label='Unfiltered Spectrum')

# # Set lims for three-level or effective two-level
# if 2 * delta == -alpha:
#     plt.xlim(-1.5 * Omega, 1.5 * Omega)
# else:
#     plt.xlim(-abs(alpha), abs(alpha))
# plt.xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=15)
# plt.ylabel(r'Power Spectrum (a.u.)', fontsize=15)

# plt.legend(loc='best', fontsize=15)
# plt.title(r'$\left( \Omega = {}, \alpha = {}, \delta = {}, \xi = {} \right), \left( \kappa = {}, \omega_{{0}} = {}, N = {}, \delta\omega = {} \right)$'.format(Omega, alpha, delta, xi, kappa, w0 ,N, dw), fontsize=15)

# plt.tight_layout()
# plt.show()

#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
# plt.figure(figsize=[10, 6])

# plt.plot(tau2, g2, color='C0', label='Filtered Correlation')
# plt.plot(tau2, g2_atom, color='k', alpha=0.5, ls='dashed', label='Unfiltered Correlation')

# plt.xlabel(r'$\gamma \tau$', fontsize=15)
# plt.ylabel(r'$g^{(2)}(\tau)$', fontsize=15)
# plt.legend(loc='best', fontsize=15)
# plt.title(r'$\left( \Omega = {}, \alpha = {}, \delta = {}, \xi = {} \right), \left( \kappa = {}, \omega_{{0}} = {}, N = {}, \delta\omega = {} \right)$'.format(Omega, alpha, delta, xi, kappa, w0 ,N, dw), fontsize=15)

# plt.tight_layout()
# plt.show()

# -----------------------------------------------------------------------------#
#                           PLOT SPECTRUM AND g^{(2)}                          #
# -----------------------------------------------------------------------------#
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 8])

# Spectrum
ax[0].plot(wlist, spec, color='C0', ls='solid', label='Filtered Spectrum')
ax[0].plot(wlist, spec_atom, color='k', ls='dashed', alpha=0.5, label='Unfiltered Spectrum')

# Set lims for three-level or effective two-level
if 2 * delta == -alpha:
    ax[0].set_xlim(-1.5 * Omega, 1.5 * Omega)
else:
    ax[0].set_xlim(-abs(alpha), abs(alpha))
ax[0].set_xlabel(r'$\left( \omega - \omega_{d} \right) / \Gamma$', fontsize=12)
ax[0].set_ylabel('Power Spectrum (a.u.)', fontsize=12)
ax[0].legend(loc='best', fontsize=12)
ax[0].set_title(r'Three-Level Atom $\left( \Omega = {} \Gamma, \alpha = {} \Gamma, \delta = {} \Gamma, \xi = {} \right)$'.format(Omega, alpha, delta, xi), fontsize=12)


# Second-order correlation
ax[1].plot(tau2, g2, color='C0', ls='solid', label='Filtered Photon Correlation')
# ax[1].plot(tau2, g2_atom, color='k', ls='dashed', alpha=0.5, label='Unfiltered Photon Correlation')

ax[1].set_xlabel(r'$\Gamma \tau$', fontsize=12)
ax[1].set_ylabel(r'$g^{(2)}(\tau)$', fontsize=12)
ax[1].set_title(r'Multi-Mode Filter $\left( N = {}, \delta\omega = {} \Gamma, \kappa = {} \Gamma, \omega_{{0}} = {} \Gamma \right)$'.format(N, dw, kappa, w0), fontsize=12)
ax[1].legend(loc='best', fontsize=12)

fig.tight_layout()
fig.show()
