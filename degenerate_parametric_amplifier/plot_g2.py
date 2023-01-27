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
kappa_p, Delta, K, N, halfwidth, kappa_f, dw, w0 = \
    np.genfromtxt("./data_files/g2_auto_parameters.txt", usecols=1, delimiter='=')
N = int(N)
    
# Read data
tau = np.genfromtxt("./data_files/g2_auto_corr.txt", usecols=0)
g2 = np.genfromtxt("./data_files/g2_auto_corr.txt", usecols=1)

# g2_full = np.genfromtxt(filename("g2_corr"), usecols=1)
# g2_unfiltered = g2_RK4(kappa_p, Delta, K, tau)

#-----------------------------------------------------------------------------#
#                                    PLOT                                     #
#-----------------------------------------------------------------------------#
fig = plt.figure(num='Auto Correlation', figsize=[8, 6])
ax = plt.gca()

if N==0:
    ax.plot(tau, g2, color='C0', ls='solid',
            label=r'Filtered: $\kappa_{{f}} = {} \kappa_{{p}}$'.format(kappa_f))
    ax.set_title((r'Single-Mode: $g^{{(2)}}(\tau)$ with $\omega_{{0}} = {} '
                  r'\kappa_{{p}}, \kappa_{{f}} = {} \kappa_{{p}}$'
                  ).format(w0, kappa_f), fontsize=12)
else:
    ax.plot(tau, g2, color='C0', ls='solid',
            label='Filtered: $\kappa_{{f}} = {} \kappa_{{p}}$'.format(kappa_f))
    ax.set_title((r'Multi-Mode $(N={}, \delta\omega = {} \kappa_{{p}})$: '
                  r'$g^{{(2)}}(\tau)$ with $\omega_{{0}} = {} '
                  r'\kappa_{{p}}, N \delta\omega = {} \kappa_{{p}}$'
                  ).format(N, dw, w0, N*dw), fontsize=12)

# ax.plot(tau, g2_unfiltered, color='C1', ls='dashed', label='Unfiltered')

ax.legend(loc='best', fontsize=12)

ax.set_xlabel(r'$\kappa_{p} \tau$', fontsize=12)
ax.set_ylabel(r'$\langle F^{\dagger} (0) F^{\dagger} F(\tau) F(0) \rangle_{ss}$', fontsize=12)

fig.tight_layout()
fig.show()