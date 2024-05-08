#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 13:25:00 2024

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../../paper_style.mplstyle')

plt.close('all')

# Parameter directorys
dir_dat = 'good_parameters'
# dir_dat = 'good_initial'

# Filenames: Right peak
dat_R_good = f"./data_files/{dir_dat}/right/g2_corr_good.txt"
par_R_good = f"./data_files/{dir_dat}/right/g2_parameters_good.txt"

dat_R_bad  = f"./data_files/{dir_dat}/right/g2_corr_bad.txt"
par_R_bad  = f"./data_files/{dir_dat}/right/g2_parameters_bad.txt"

# Filenames: Central peak
dat_C_good = f"./data_files/{dir_dat}/centre/g2_corr_good.txt"
par_C_good = f"./data_files/{dir_dat}/centre/g2_parameters_good.txt"

dat_C_bad  = f"./data_files/{dir_dat}/centre/g2_corr_bad.txt"
par_C_bad  = f"./data_files/{dir_dat}/centre/g2_parameters_bad.txt"

# Filename: Figure
filename_out = "../../../images/sect5/fig11_combined_g2_good_vs_bad.pdf"

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def g2_dressed_states(tau_in, Omega_in, gamma_in, w0_in):
    """
    Analytic expression for the normalised second-order correlation function
    of the dressed state transitions.
    """
    from numpy import exp, ones
    if abs(w0_in) == abs(Omega_in):
        g2_f = 1.0 - exp(-0.5 * gamma_in * tau_in)
    elif abs(w0_in) == 0.0:
        g2_f = ones(len(tau_in))
    return g2_f

#-----------------------------------------------------------------------------#
#                                  DATA THINGS                                #
#-----------------------------------------------------------------------------#
# Parameters
gamma, Omega, N, halfwidth, kappa, dw, m, w0 = \
    np.genfromtxt(par_R_good, delimiter="=", usecols=1)
N = int(N)    
del halfwidth, kappa, dw

# Read tau
tau = np.genfromtxt(dat_R_good, usecols=0)

# Dressed state correlation
g2_R = g2_dressed_states(tau, Omega, gamma, Omega)
g2_C = g2_dressed_states(tau, Omega, gamma, 0.0)

#-------------------------------#
#     Read Data: Right Peak     #
#-------------------------------#
# Single-mode data
g2_R_good_single = np.genfromtxt(dat_R_good, usecols=1)
g2_R_bad_single  = np.genfromtxt(dat_R_bad, usecols=1)

# Multi-mode data
g2_R_good_multi  = np.genfromtxt(dat_R_good, usecols=2)
g2_R_bad_multi   = np.genfromtxt(dat_R_bad, usecols=2)

# Parameters: Bad
K_R_bad, kappa_R_bad, dw_R_bad = \
  np.genfromtxt(par_R_bad, delimiter='=', usecols=1, skip_header=3, skip_footer=2)

# Parameters: Good
K_R_good, kappa_R_good, dw_R_good = \
  np.genfromtxt(par_R_good, delimiter='=', usecols=1, skip_header=3, skip_footer=2)

#---------------------------------#
#     Read Data: Central Peak     #
#---------------------------------#
# Single-mode data
g2_C_good_single = np.genfromtxt(dat_C_good, usecols=1)
g2_C_bad_single  = np.genfromtxt(dat_C_bad, usecols=1)

# Multi-mode data
g2_C_good_multi  = np.genfromtxt(dat_C_good, usecols=2)
g2_C_bad_multi   = np.genfromtxt(dat_C_bad, usecols=2)

# Parameters: Bad
K_C_bad, kappa_C_bad, dw_C_bad = \
  np.genfromtxt(par_C_bad, delimiter='=', usecols=1, skip_header=3, skip_footer=2)

# Parameters: Good
K_C_good, kappa_C_good, dw_C_good = \
  np.genfromtxt(par_C_good, delimiter='=', usecols=1, skip_header=3, skip_footer=2)

#--------------------------#
#     Print Halfwidths     #
#--------------------------#
print("(Right) Good Halfwidth : {}".format(K_R_good))
print("(Right) Bad Halfwidth  : {}".format(K_R_bad))

print("(Centre) Good Halfwidth: {}".format(K_C_good))
print("(Centre) Bad Halfwidth : {}".format(K_C_bad))

#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
fig, ax = plt.subplots(num='Good vs Bad', sharex=False, nrows=2, ncols=1,
                       figsize=[6, 6])

#----------------------#
#     Central Peak     #
#----------------------#
# Plot multi-mode
ax[0].plot(tau, g2_C_good_multi, color='C0', ls='solid',
          label=r'Multi-Mode $\left( K = {:.1f} \gamma \right)$'.format(K_C_good))
# Plot single-mode
ax[0].plot(tau, g2_C_bad_single, color='C1', ls='dashed', alpha=1.0,
           label=r'Single-Mode $\left( K = {:.1f} \gamma \right)$'.format(K_C_bad))
# Plot dressed state approximation
ax[0].plot(tau, g2_C, color='k', alpha=0.5, ls='dotted',
          label='Eq.~(38a)')

#--------------------#
#     Right Peak     #
#--------------------#
# Plot multi-mode
ax[1].plot(tau, g2_R_good_multi, color='C0', ls='solid',
          label=r'Multi-Mode $\left( K = {:.1f} \gamma \right)$'.format(K_R_good))
# Plot single-mode
ax[1].plot(tau, g2_R_bad_single, color='C1', ls='dashed', alpha=1.0,
           label=r'Single-Mode $\left( K = {:.1f} \gamma \right)$'.format(K_R_bad))
# Plot dressed state approximation
ax[1].plot(tau, g2_R, color='k', alpha=0.5, ls='dotted',
           label='Eq.~(38c)')

#------------------#
#     Add Text     #
#------------------#
# Central peak
ax[0].text(x=-0.7, y=1.15, s='(a)')
ax[1].text(x=-0.7, y=1.0, s='(b)')

#--------------------#
#     Axis Stuff     #
#--------------------#
# Cycle through axes
for j in range(2):
    # Legend
    ax[j].legend(loc='lower right')
    
    # Grid
    ax[j].grid()
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax[j].set_xticks(np.arange(0.0, 6, 1))
    ax[j].set_xticks(np.arange(0.5, 6, 1), minor=True)

    if j == 0:
        ax[j].set_yticks(np.arange(0.90, 1.20, 0.05))
        ax[j].set_yticks(np.arange(0.775, 1.5, 0.15), minor=True)
    elif j == 1:
        ax[j].set_yticks(np.arange(0.0, 1.25, 0.25))
        ax[j].set_yticks(np.arange(0.125, 1.25, 0.25), minor=True)

    # ax[i, j].set_yticks(np.arange(0.0, 1.2, 0.2))
    # ax[i, j].set_yticks(np.arange(0.1, 1.1, 0.2), minor=True)

    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax[j].set_xlim(-0.05, 5.05)

    if j == 0:
        ax[j].set_ylim(0.89, 1.16)
    elif j == 1:
        ax[j].set_ylim(-0.05, 1.05)

    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax[j].set_xlabel(r'$\gamma \tau$')
    if j == 0:
        ax[j].set_ylabel(r'$g^{(2)}(0, 0; 0, \tau)$')
    elif j == 1:
        ax[j].set_ylabel(r'$g^{(2)}(\Omega, 0; \Omega, \tau)$')

#----------------------#
#     Figure Stuff     #
#----------------------#
# fig.tight_layout(pad=0)
fig.savefig(filename_out)
fig.show()