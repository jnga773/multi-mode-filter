# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

plt.close('all')

def filename(file_in):
    direct_f = "./data_files/"
    ext_f = ".txt"
    filename_out = direct_f + file_in + ext_f
    return filename_out

#------------------------------------------------------------------------------#
#                                FILENAME THINGS                               #
#------------------------------------------------------------------------------#
# Read parameters
gamma, Omega, kappa, dw, epsilon, N, phase = \
    np.genfromtxt(filename("parameters"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
Omega = round(Omega, 2)

# Pull data from file
# Central frequencies
w0_list = np.genfromtxt(filename("g2_0"), usecols=0)
# Steady state photon number
g2_0 = np.genfromtxt(filename("g2_0"), usecols=1)

#------------------------------------------------------------------------------#
#                            PLOT PHOTON NUMBER SCAN                           #
#------------------------------------------------------------------------------#
plt.figure(figsize=[12, 8])

# Plot photon number scan
plt.plot(w0_list, g2_0)

plt.xlim(-1.5 * Omega, 1.5 * Omega)

plt.ylabel(r'$ \langle A^{\dagger} A^{\dagger} A A \rangle_{ss} / \langle A^{\dagger} A \rangle^{2}_{ss} $', fontsize=15)
# plt.ylabel(r'$ \langle A^{\dagger} A^{\dagger} A A \rangle_{ss} $', fontsize=15)
plt.xlabel(r'$ \omega_{0} / \gamma $', fontsize=15)

if N > 0:
    plt.title(r'$g^{{(2)}}(\tau=0)$ with $\left( \Omega = {}, \kappa = {} \gamma, \delta\omega = {} \gamma, N = {} \right)$'.format(Omega, kappa, dw, N), fontsize=12)
else:
    plt.title(r'$g^{{(2)}}(\tau=0)$ with $\left( \Omega = {} \gamma, \kappa = {} \gamma \right)$'.format(Omega, kappa), fontsize=12)
    # plt.title(r'$\langle a^{{\dagger}} a^{{\dagger}} a a \rangle_{{ss}}$ with $\left( \Omega = {} \gamma, \kappa = {} \gamma \right)$'.format(Omega, kappa), fontsize=12)
   
plt.tight_layout()
plt.show()
