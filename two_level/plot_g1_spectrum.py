# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
# from scipy.signal import find_peaks
import sys
# Append path to add Jacobs_Function
sys.path.append("/home/jnga773/Documents/898-PhD/test/")
# sys.path.append("E:/Google Drive/University/Physics/898 PhD/test")
from Jacobs_Functions import filename, spectrum

plt.close('all')

def mollow_triplet(tau_in, gamma_in, Omega_in):
    """
    Expression from Howard vol. 1
    """
    # Turn inputs into complex values
    Omega_in = np.complex(Omega_in)
    gamma_in = np.complex(gamma_in)
    # input tau
    # Calculate first order correlation
    Yf = (np.sqrt(2.0) * Omega_in) / (gamma_in)
    df = np.sqrt(((0.25 * gamma_in) ** 2) - (Omega_in ** 2))
    g1_f = (0.25 * ((Yf ** 2) / (1 + (Yf ** 2)))) * np.exp(-0.5 * gamma_in * tau_in) - \
           (0.125 * ((Yf ** 2) / ((1 + (Yf ** 2)) ** 2))) * (1 - (Yf ** 2) + ((1 - 5 * (Yf ** 2)) * (0.25 * gamma_in) / df)) * np.exp(-((0.75 * gamma_in) - df) * tau_in) - \
           (0.125 * ((Yf ** 2) / ((1 + (Yf ** 2)) ** 2))) * (1 - (Yf ** 2) - ((1 - 5 * (Yf ** 2)) * (0.25 * gamma_in) / df)) * np.exp(-((0.75 * gamma_in) + df) * tau_in)
    print("Mean Im[g1_atom] = {}".format(np.mean(np.imag(g1_f))))
    return g1_f

# Normalisation function
def norm_spectra(Omega_in, atom_spec_in, filtered_spec_in, w0_in):
    """
    Depending on where the filter is centred, normalise the filtered spectrum
    to the relevant peak.
    """
    from numpy import array
    # Calculate where the peaks are located in the UNfiltered spectrum
    from scipy.signal import find_peaks
    peaks = find_peaks(atom_spec_in, height=0.1)[1]["peak_heights"]
    
    # Calculate where the peaks should fall from eigenvalues
    peaks_calculated = array([-Omega_in, 0.0, Omega_in])
    # print(peaks_calculated)
    
    # Cycle through peaks_calculated and compare with w0
    norm_peak_value = 1.0
    for place, item in enumerate(peaks_calculated):
        if round(item, 2) == round(w0_in, 2):
            # save the peak value
            norm_peak_value = peaks[place]
            # Leave loop
            break
    # Renormalise spectra
    spec_out = filtered_spec_in * norm_peak_value
    return spec_out

#------------------------------------------------------------------------------#
#                                FILENAME THINGS                               #
#------------------------------------------------------------------------------#
# Read parameters
gamma, Omega, w0, kappa, dw, epsilon, N, phase, dt, tau1_max = \
    np.genfromtxt(filename("parameters"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
Omega = round(Omega, 2)
w0 = round(w0, 2)

# Pull data from file
tau = np.genfromtxt(filename("g1_corr"), usecols=0)
corr = np.genfromtxt(filename("g1_corr"), usecols=1) + 1j * \
        np.genfromtxt(filename("g1_corr"), usecols=2)
# tau = np.genfromtxt(filename("me_g1_corr"), usecols=0)
corr_atom = mollow_triplet(tau, gamma, Omega)

# Calculate Fourier transform
spec, wlist = spectrum(tau, corr, norm='peak')
spec_atom, wlist = spectrum(tau, corr_atom, norm='peak')

# Renormalise filtered spectra
spec = norm_spectra(Omega, spec_atom, spec, w0)

#-----------------------------------------------------------------------------#
#                               PLOT G^{(1)}                                  #
#-----------------------------------------------------------------------------#
# fig, ax = plt.subplots(2, 1, sharex=True, figsize=[8, 6])

# # Real part
# ax[0].plot(tau, corr.real, color='C0', ls='solid', lw=2.0, label='Filtered Spectrum')
# ax[0].plot(tau, corr_atom.real, color='C1', ls='dashed', lw=1.0, label='Full Spectrum')
# ax[0].set_ylabel(r'Real Part', fontsize=15)
# ax[0].set_xlim(-0.2, 5.2)
# ax[0].legend(loc='best', fontsize=15)
# ax[0].set_title(r'$G^{{(1)}}(\tau)$ with $\left( \Omega = {} \gamma \right)$'.format(Omega), fontsize=15)

# ax[1].plot(tau, corr.imag, color='C0', ls='solid', lw=2.0)
# ax[1].plot(tau, corr_atom.imag, color='k', ls='dashed', lw=1.0, alpha=0.5)
# ax[1].set_ylabel(r'Imaginary Part', fontsize=15)
# ax[1].set_xlabel(r'$\gamma t$', fontsize=15)
# ax[1].set_title(r'$N = {}, \kappa = {} \gamma, \omega_{{0}} = {} \gamma, \delta\omega = {} \gamma, \epsilon = {}$'.format(N, kappa, w0, dw, epsilon), fontsize=15)


# fig.tight_layout()
# fig.show()

#-----------------------------------------------------------------------------#
#                               PLOT SPECTRUM                                 #
#-----------------------------------------------------------------------------#
# plt.figure(figsize=[8, 8])

# plt.plot(wlist, spec, color='C0', ls='solid', lw=2.0, label='Filtered Spectrum')
# plt.plot(wlist, spec_atom, color='k', ls='dashed', lw=1.0, alpha=0.5, label='Full Spectrum')

# plt.xlim(-1.2 * Omega, 1.2 * Omega)
# plt.xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=15)
# plt.ylabel(r'Power Spectrum (a.u.)', fontsize=15)

# plt.legend(loc='best', fontsize=15)
# plt.title(r'$\Omega = {} \gamma, N = {}, \delta\omega = {} \gamma, \kappa = {} \gamma, \omega_{{0}} = {} \gamma$'.format(Omega, N, dw, kappa, w0))

# plt.tight_layout()
# plt.show()

#-----------------------------------------------------------------------------#
#                          PLOT G^{(1)} AND SPECTRUM                          #
#-----------------------------------------------------------------------------#
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 8])

# First-order correlation
ax[0].plot(tau, corr.real, color='C0', ls='solid', lw=2.0, label='Real')
ax[0].plot(tau, corr.imag, color='C1', ls='dashed', lw=1.0, label='Imaginary')
ax[0].set_xlim(-0.1, 10.1)

ax[0].set_xlabel(r'$\gamma \tau$', fontsize=12)
ax[0].set_ylabel(r'$G^{(1)}(\tau)$', fontsize=12)
ax[0].set_title(r'$\Omega = {} \gamma$'.format(Omega), fontsize=12)
ax[0].legend(loc='best', fontsize=12)

# Spectra
ax[1].plot(wlist, spec, color='C0', ls='solid', lw=2.0, label='Filtered')
ax[1].plot(wlist, spec_atom, color='k', ls='dashed', lw=1.0, alpha=0.5, label='Unfiltered')

ax[1].set_xlim(-1.2 * Omega, 1.2 * Omega)
ax[1].set_xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=12)
ax[1].set_ylabel('Power Spectrum (a.u.)', fontsize=12)
ax[1].set_title(r'$N = {}, \delta\omega = {} \gamma, \kappa = {} \gamma, \omega_{{0}} = {} \gamma$'.format(N, dw, kappa, w0), fontsize=12)
ax[1].legend(loc='best', fontsize=12)

fig.tight_layout()
fig.show()