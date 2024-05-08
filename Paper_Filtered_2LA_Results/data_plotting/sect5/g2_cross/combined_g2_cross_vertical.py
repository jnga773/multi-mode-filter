#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 12:54:53 2024

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../../paper_style.mplstyle')

plt.close('all')

# Directories: Right to central peak
dir_RC = './data_files/centre/'
# Directories: Right to left peak
dir_RL = './data_files/right/'

# Figure filename
filename_out = '../../../images/sect5/fig14_combined_g2_cross.pdf'

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def calc_dressed_g2(gamma_in, Omega_in, tau_in, halfwidth_in, which_peak_in):
    if which_peak_in == "centre":
        # g2_out = np.ones(len(tau_in))
        
        g2_out = (1 - np.exp(-halfwidth_in * tau_in)) ** 2
        
    elif which_peak_in == "right":
        # g2_out = 1 + np.exp(-0.5 * gamma_in * tau_in)
        
        g2_out = (np.exp(-0.5 * gamma_in * tau_in) - 1)
        g2_out += 2 * ((1 - 0.5*np.exp(-halfwidth_in * tau_in)) ** 2 )
        g2_out += 0.5 * np.exp(-2 * halfwidth_in * tau_in)
        
    return g2_out

#-----------------------------------------------------------------------------#
#                                  DATA THINGS                                #
#-----------------------------------------------------------------------------#
#------------------------------------------#
#     Read Data: Right to Central Peak     #
#------------------------------------------#
# Read parameters
gamma, Omega, N, K_RC, kappa_RC, dw = \
    np.genfromtxt(dir_RC + 'g2_parameters.txt', delimiter="=", usecols=1, max_rows=6)
N = int(N)

# Tau time
tau          = np.genfromtxt(dir_RC + 'g2_corr.txt', usecols=0)
# Single-mode
g2_RC_single = np.genfromtxt(dir_RC + 'g2_corr.txt', usecols=1)
# Multi-mode
g2_RC_multi  = np.genfromtxt(dir_RC + 'g2_corr.txt', usecols=2)

# Calculate secular approximation
g2_RC_multi_dressed  = calc_dressed_g2(gamma, Omega, tau, K_RC, 'centre')
g2_RC_single_dressed = calc_dressed_g2(gamma, Omega, tau, kappa_RC, 'centre')

#------------------------------------------#
#     Read Data: Right to Central Peak     #
#------------------------------------------#
# Read parameters
gamma, Omega, N, K_RL, kappa_RL, dw = \
    np.genfromtxt(dir_RL + 'g2_parameters.txt', delimiter="=", usecols=1, max_rows=6)
N = int(N)
# Tau time
# tau          = np.genfromtxt(dir_RC + 'g2_corr.txt', usecols=0)
# Single-mode
g2_RL_single = np.genfromtxt(dir_RL + 'g2_corr.txt', usecols=1)
# Multi-mode
g2_RL_multi  = np.genfromtxt(dir_RL + 'g2_corr.txt', usecols=2)

# Calculate secular approximation
g2_RL_multi_dressed  = calc_dressed_g2(gamma, Omega, tau, K_RL, 'right')
g2_RL_single_dressed = calc_dressed_g2(gamma, Omega, tau, kappa_RL, 'right')

#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
fig, ax = plt.subplots(num="g2 cross (Centre and Right)", nrows=2, ncols=1,
                       sharex=False, figsize=[6, 6])

#----------------------------#
#     Plot: Central Peak     #
#----------------------------#
# Plot: Multi-Mode
ax[0].plot(tau, g2_RC_multi, color='C0', ls='solid',
           label=rf'Multi-Mode $\left( K = {K_RC:.2f} \gamma \right)$')
# Plot: Single-mode
ax[0].plot(tau, g2_RC_single, color='C1', ls='dashed',
           label=rf'Single-Mode $\left( K = {kappa_RC:.2f} \gamma \right)$')

# Plot: Secular approximations
ax[0].plot(tau, g2_RC_multi_dressed, color='k', ls='dotted', alpha=0.5,
           label='Eq.~(40a)')
ax[0].plot(tau, g2_RC_single_dressed, color='k', ls='dotted', alpha=0.5)

# Legend
ax[0].legend(loc='lower right')

#--------------------------#
#     Plot: Right-Peak     #
#--------------------------#
# Plot: Multi-Mode
ax[1].plot(tau, g2_RL_multi, color='C0', ls='solid',
           label=rf'Multi-Mode $\left( K = {K_RL:.2f} \gamma \right)$')
# Plot: Single-mode
ax[1].plot(tau, g2_RL_single, color='C1', ls='dashed',
           label=rf'Single-Mode $\left( K = {kappa_RL:.2f} \gamma \right)$')

# Plot: Secular approximations
ax[1].plot(tau, g2_RL_multi_dressed, color='k', ls='dotted', alpha=0.5,
           label='Eq.~(40b)')
ax[1].plot(tau, g2_RL_single_dressed, color='k', ls='dotted', alpha=0.5)

# Legend
ax[1].legend()

#-------------------------#
#     Add Figure Text     #
#-------------------------#
ax[0].text(x=-0.65, y=1.25, s='(a)')
ax[1].text(x=-0.65, y=2.00, s='(b)')

#--------------------#
#     Axis Ticks     #
#--------------------#
# Ticks: X-axis
ax[0].set_xticks(np.arange(0.0, 6.0, 1.0))
ax[0].set_xticks(np.arange(0.5, 6.0, 1.0), minor=True)

ax[1].set_xticks(np.arange(0.0, 6.0, 1.0))
ax[1].set_xticks(np.arange(0.5, 6.0, 1.0), minor=True)

# Ticks: Y-axis
ax[0].set_yticks(np.arange(0.0, 2.0, 0.25))
ax[0].set_yticks(np.arange(0.125, 2.0, 0.25), minor=True)

ax[1].set_yticks(np.arange(0.0, 2.25, 0.25))
ax[1].set_yticks(np.arange(0.125, 2.25, 0.25), minor=True)

#---------------------#
#     Axis Limits     #
#---------------------#
# X-axis
ax[0].set_xlim(-0.05, 5.05)
ax[1].set_xlim(-0.05, 5.05)

# Y-axis
ax[0].set_ylim(-0.02, 1.27)
ax[1].set_ylim(0.93, 2.02)

#---------------#
#     Grids     #
#---------------#
ax[0].grid()
ax[1].grid()


#---------------------#
#     Axis Labels     #
#---------------------#
ax[0].set_xlabel(r'$\gamma \tau$')
ax[0].set_ylabel(r'$g^{(2)}(\Omega, 0; 0, \tau)$')

ax[1].set_xlabel(r'$\gamma \tau$')
ax[1].set_ylabel(r'$g^{(2)}(\Omega, 0; -\Omega, \tau)$')

#----------------------#
#     Figure Stuff     #
#----------------------#
fig.tight_layout()
fig.savefig(filename_out)
fig.show()
