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

def filename(filename_in, parent_in="./data_files/", ext_in="txt"):
    filename_out = parent_in + filename_in + "." + ext_in
    return filename_out

#------------------------------------------------------------------------------#
#                                FILENAME THINGS                               #
#------------------------------------------------------------------------------#
# Read parameters
gamma, Omega, w0, kappa, dw, epsilon, N, phase, dt, t_max = \
    np.genfromtxt(filename("parameters"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
Omega = round(Omega, 2)
w0 = round(w0, 2)

# Read data
t = np.genfromtxt(filename("states"), usecols=0)
gg = np.genfromtxt(filename("states"), usecols=1)
ee = np.genfromtxt(filename("states"), usecols=2)
photon = np.genfromtxt(filename("states"), usecols=3)
sz = ee - gg
# Tried and tested master equation
gg_me = np.genfromtxt(filename("old_states"), usecols=1)
ee_me = np.genfromtxt(filename("old_states"), usecols=2)
photon_me = np.genfromtxt(filename("old_states"), usecols=3)
sz_me = ee_me - gg_me

#-----------------------------------------------------------------------------#
#                               PLOT STATES                                   #
#-----------------------------------------------------------------------------#

# fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=[8, 6])

# # < \sigma_{z}(t) >
# ax[0].plot(t, sz, color='C0', ls='solid', lw=2.0, label='Newer')
# ax[0].plot(t, sz_me, color='C1', ls='dashed', lw=1.0, label='Older')
# ax[0].set_ylim(-1.05, 1.05)
# ax[0].set_ylabel(r'$\langle \hat{\sigma}_{z}(t) \rangle$', fontsize=15)
# ax[0].set_title(r'Two-Level Atom $\left( \Omega = {} \gamma \right)$'.format(round(Omega, 2)), fontsize=15)

# # < A^{\dagger} A(t) >
# ax[1].plot(t, photon, color='C0', ls='solid', lw=2.0)
# ax[1].plot(t, photon_me, color='C1', ls='dashed', lw=1.0)
# ax[1].set_ylabel(r'$\langle \hat{A}^{\dagger} \hat{A}(t) \rangle$', fontsize=15)
# ax[1].set_title(r'$N = {}, \kappa = {} \gamma, \omega_{{0}} = {} \gamma, \delta\omega = {} \gamma, \epsilon = {}$'.format(N, kappa, w0, dw, epsilon), fontsize=15)

# # Labels
# ax[0].legend(loc='best', fontsize=15)
# ax[1].set_xlabel(r'$\gamma t$', fontsize=15)

# fig.tight_layout()
# fig.show()

#-----------------------------------------------------------------------------#
#                              PLOT MOMENTS                                   #
#-----------------------------------------------------------------------------#
# Operator moments
moment = np.genfromtxt(filename("operators"), usecols=1) + 1j * \
          np.genfromtxt(filename("operators"), usecols=2)
moment_me = np.genfromtxt(filename("old_operators"), usecols=1) + 1j * \
            np.genfromtxt(filename("old_operators"), usecols=2)
            
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=[8, 6])

# Real part
ax[0].plot(t, moment.real, color='C0', ls='solid', lw=2.0, label='Newer')
ax[0].plot(t, moment_me.real, color='C1', ls='dashed', lw=1.0, label='Older')
ax[0].set_ylabel(r'Real Part', fontsize=15)
ax[0].legend(loc='best', fontsize=15)
ax[0].set_title(r'$N = {}, \kappa = {} \gamma, \omega_{{0}} = {} \gamma, \delta\omega = {} \gamma, \epsilon = {}$'.format(N, kappa, w0, dw, epsilon), fontsize=15)


ax[1].plot(t, moment.imag, color='C0', ls='solid', lw=2.0)
ax[1].plot(t, moment_me.imag, color='C1', ls='dashed', lw=1.0)
ax[1].set_ylabel(r'Imaginary Part', fontsize=15)
ax[1].set_xlabel(r'$\gamma t$', fontsize=15)

fig.tight_layout()
fig.show()
