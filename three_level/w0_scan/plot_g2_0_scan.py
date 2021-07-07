# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

from _my_functions import filename

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

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
w0_list = np.genfromtxt(filename("g2_0"), usecols=0)
# Steady state photon number
g2_0 = np.genfromtxt(filename("g2_0"), usecols=1)

#------------------------------------------------------------------------------#
#                            PLOT PHOTON NUMBER SCAN                           #
#------------------------------------------------------------------------------#
plt.figure(figsize=[12, 8])

# Plot photon number scan
plt.plot(w0_list, g2_0)

plt.ylabel(r'$ g^{(2)}(\tau=0) $', fontsize=15)
plt.xlabel(r'$ \omega_{0} / \gamma $', fontsize=15)
# plt.title(r'$ \Omega = %s, \mathrm{phase} = %s\pi, \kappa = %s, \delta\omega = %s, N = %s $'%(Omega, phase, kappa, dw, N), fontsize=15)

plt.tight_layout()
plt.show()
