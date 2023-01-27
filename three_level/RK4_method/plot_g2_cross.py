# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

from Jacobs_Functions import filename
from analytic_functions._dressed_state_correlations import dressed_g2_calc, three_level_eig

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

plt.close('all')

#-----------------------------------------------------------------------------#
#                               FILENAME THINGS                               #
#-----------------------------------------------------------------------------#
# Read parameters
Gamma, Omega, alpha, delta, xi, N, halfwidth, kappa, dw, w0a, w0b = \
    np.genfromtxt(filename("g2_cross_parameters"), delimiter="=", usecols=1)
N = int(N)

# Pull data from file
tau = np.genfromtxt(filename("g2_cross_multi"), usecols=0)
g2 = np.genfromtxt(filename("g2_cross_multi"), usecols=1)

g2_single = np.genfromtxt(filename("g2_cross_single"), usecols=1)

dt = 0.001
tau_half = np.arange(0, tau[-1] + dt, dt)

# Half data
# tau_half = np.genfromtxt(filename("positive"), usecols=0)
# pos = np.genfromtxt(filename("positive"), usecols=1)
# neg = np.genfromtxt(filename("negative"), usecols=1)

# g2 = np.concatenate((np.flip(neg), pos[1:]))
# tau = np.concatenate((-np.flip(tau_half), tau_half[1:]))

dressed_neg = dressed_g2_calc(tau_half, w0b, Gamma, Omega, alpha, delta, xi, w0b_in=w0a)
dressed_pos = dressed_g2_calc(tau_half, w0a, Gamma, Omega, alpha, delta, xi, w0b_in=w0b)
g2_dressed = np.concatenate((np.flip(dressed_neg), dressed_pos[1:]))

#-----------------------------------------------------------------------------#
#                           PLOT CROSS-CORRELATION                            #
#-----------------------------------------------------------------------------#
fig = plt.figure(figsize=[10, 6])
ax = plt.gca()

# ax.plot(tau, g2, color='C0',
#         label=(r'Filtered Correlation $\left( N = {}, \delta\omega = {}'
#                r'\Gamma, \kappa = {} \Gamma \right)$').format(N, round(dw, 2), kappa))
ax.plot(tau, g2, color='C0', ls='solid',
        label='Improved Filtering Method (halfwidth = $ 8 \Gamma $)')

ax.plot(tau, g2_single, color='C1', ls='dashed',
        label='Single-Mode Filter (halfwidth = $ 4 \Gamma $)')

# Plot dressed-state correlation
ax.plot(tau, g2_dressed, color='k', ls='dashed', alpha=0.5, 
        label='Dressed State Correlation')

# Move the y-axis to the centre
ax.spines['left'].set_position('zero')
ax.spines.top.set_visible(False)
ax.spines.right.set_visible(False)

# Change tick direction
ax.tick_params(axis='y', direction='in')

ax.set_xlabel(r'$\Gamma \tau$', fontsize=12)
# ax.set_ylabel(r'$g^{(2)}_{\mathrm{Auto}}(\tau)$', fontsize=15)

# ax.set_title((r'$g^{{(2)}}_{{\mathrm{{cross}}}}(\tau)$ with $\left( \Omega = {}'
#               r'\Gamma, \alpha = {}, \delta = {}, \xi = {}'
#               r'\right)$ and $\left( \omega_{{0}}^{{(A)}} = {} \Gamma,'
#               r' \omega_{{0}}^{{(B)}} = {} \Gamma \right)$'
#               ).format(round(Omega, 2), alpha, delta, xi, round(w0a, 2), round(w0b, 2)),
#              fontsize=12)

ax.legend(loc='best', fontsize=12)

fig.tight_layout()
# fig.savefig("./koala_abstract_2022.pdf")
fig.savefig("./koala_abstract_2022.png", dpi=1000)
fig.show()