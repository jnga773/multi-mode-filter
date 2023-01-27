#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 14:24:06 2022

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

# from parametric_functions import g2_RK4

plt.close('all')

#-----------------------------------------------------------------------------#
#                               FILENAME THINGS                               #
#-----------------------------------------------------------------------------#
# Grab parameters
kappa_p, Delta, K, N, halfwidth, kappa_f, dw, w0a, w0b = \
    np.genfromtxt("./data_files/g2_cross_parameters.txt", usecols=1, delimiter='=')
N = int(N)
    
# Read data
tau = np.genfromtxt("./data_files/g2_cross_corr.txt", usecols=0)

# Auto-correlation
g2_auto = np.genfromtxt("./data_files/g2_auto_corr.txt", usecols=1)

# Cross-correlation
g2_cross = np.genfromtxt("./data_files/g2_cross_corr.txt", usecols=1)

print("G2(0) (auto) = {}".format(g2_auto[0]))
print("G2(0) (cross) = {}".format(g2_cross[0]))
# print(' ')
# print("G2(-1) (OLD) = {}".format(g2_old[-1]))
# print("G2(-1) (NEW) = {}".format(g2_sub[-1]))

#-----------------------------------------------------------------------------#
#                                    PLOT                                     #
#-----------------------------------------------------------------------------#
fig = plt.figure(num="cross vs auto", figsize=[8, 6])
ax = plt.gca()

# Plot single-mode
ax.plot(tau, g2_auto, color='C0', ls='solid',
        label=r'Single-mode $( \kappa_{{f}} = {} )$'.format(kappa_f))

# Plot mulit-mode
ax.plot(tau, g2_cross, color='C1', ls='dashed',
        label=r'Multi-mode $(N = {}, \delta\omega = {}, \kappa_{{f}} = {} )$'.format(N, dw, kappa_f))

# Labels
ax.legend(loc='best', fontsize=12)
ax.set_xlabel(r'$\kappa_{p} \tau$', fontsize=12)
ax.set_ylabel(r'$g^{(2)}(\tau)$', fontsize=12)

ax.set_title((r'Auto-Correlation with $\left( \lambda = {}, \omega_{{0}}^{{(f)}}'
              r' = \omega_{{0}}^{{(g)}} = {} \right)$').format(K, w0a))

fig.tight_layout()
fig.show()