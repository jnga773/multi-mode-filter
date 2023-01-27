#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:05:09 2022

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
g2_cross = np.genfromtxt("./data_files/g2_cross_corr.txt", usecols=1)

# g2_single = np.genfromtxt(filename("g2_cross_single"), usecols=1)

#-----------------------------------------------------------------------------#
#                                    PLOT                                     #
#-----------------------------------------------------------------------------#
fig = plt.figure(num='Cross Correlation', figsize=[6, 4])
ax = plt.gca()  

if N > 0:
    ax.plot(tau, g2_cross, color='C0', ls='solid',
            label=(r'Multi-Mode: $N = {}, \delta\omega = {}, \kappa_{{f}} = {}$'
                    ).format(N, dw, kappa_f))
else:
    ax.plot(tau, g2_cross, color='C0', ls='solid',
            label=r'Single-Mode: $\kappa_{{f}} = {} \kappa_{{p}}$'.format(kappa_f))

ax.legend(loc='best', fontsize=12)

ax.set_xlabel(r'$\kappa_{p} \tau$', fontsize=12)
ax.set_ylabel((r'$\langle F^{\dagger} (0) G^{\dagger} G(\tau) F(0) \rangle_{ss}'
               r' / \langle F^{\dagger} F \rangle_{ss} \langle G^{\dagger} '
               r'G \rangle_{ss} $'
               ), fontsize=12)

ax.set_title((r'Cross-Correlation with $\left( \lambda = {} \kappa_{{p}}, '
              r'\omega_{{0}}^{{(f)}} = {} \kappa_{{p}}, \omega_{{0}}^{{(g)}} '
              r'= {} \kappa_{{p}} \right)$'
              ).format(K, w0a, w0b), fontsize=12)

fig.tight_layout()
fig.show()