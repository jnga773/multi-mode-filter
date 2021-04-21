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
from Jacobs_Functions import filename, spectrum

plt.close('all')

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

# Pull data from file
# Central frequencies
w0_list = np.genfromtxt(filename("photon"), usecols=0)
# Steady state photon number
photon = np.genfromtxt(filename("photon"), usecols=1)

# Atomic/unfiltered correlation
tau = np.genfromtxt(filename("atom_g1_corr"), usecols=0)
g1_atom = np.genfromtxt(filename("atom_g1_corr"), usecols=1) + 1j * \
          np.genfromtxt(filename("atom_g1_corr"), usecols=2)

# Calculate Fourier transform
spec_atom, wlist = spectrum(tau, g1_atom, norm='peak')

# Renormalise atomic spectra
spec_atom = max(photon) * spec_atom / max(spec_atom)

#------------------------------------------------------------------------------#
#                            PLOT PHOTON NUMBER SCAN                           #
#------------------------------------------------------------------------------#
plt.figure(figsize=[12, 8])

# Plot photon number scan
plt.plot(w0_list, photon)
# Plot bare spectrum of atom
plt.plot(wlist, spec_atom, color='k', ls='dotted', alpha=0.5)

plt.xlim(w0_list.min(), w0_list.max())

plt.ylabel(r'$ \langle A^{\dagger} A \rangle_{ss} $', fontsize=15)
plt.xlabel(r'$ \omega_{0} / \gamma $', fontsize=15)
# plt.title(r'$ \Omega = %s, \mathrm{phase} = %s\pi, \kappa = %s, \delta\omega = %s, N = %s $'%(Omega, phase, kappa, dw, N), fontsize=15)

plt.tight_layout()
plt.show()
