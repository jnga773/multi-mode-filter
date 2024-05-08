#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 15:49:02 2024

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../../paper_style.mplstyle')

plt.close('all')

#----------------------#
#     Central Peak     #
#----------------------#
# Parameters: Multi-mode
par_C_multi  = './data_files/centre/g2_initial_parameters_multi.txt'
# # Parameters: Single-mode
# par_C_single = './data_files/centre/g2_initial_parameters_single.txt'

# Data: Multi-mode
dat_C_multi  = './data_files/centre/g2_initial_multi.txt'
# Data: Single-mode
dat_C_single = './data_files/centre/g2_initial_single.txt'

#--------------------#
#     Right Peak     #
#--------------------#
# # Parameters: Multi-mode
# par_R_multi  = './data_files/right/g2_initial_parameters_multi.txt'
# # Parameters: Single-mode
# par_R_single = './data_files/right/g2_initial_parameters_single.txt'

# Data: Multi-mode
dat_R_multi  = './data_files/right/g2_initial_multi.txt'
# Data: Single-mode
dat_R_single = './data_files/right/g2_initial_single.txt'

#-------------------------#
#     Figure Filename     #
#-------------------------#
filename_out = '../../../images/sect5/fig12_g2_initial_multi_single.pdf'

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def calc_dressed_g2(gamma_in, Omega_in, w0_in, tau_in):
    if w0_in == 0.0:
        g2_out = np.ones(len(tau_in))
    elif w0_in == Omega_in or w0_in == -Omega_in:
        g2_out = 1 - np.exp(-0.5 * gamma_in * tau_in)
        
    return g2_out

def calc_min(g2_initial_filtered_in, widths_in, min_or_max='min'):
    """
    Calculates the minimum value of g2_initial_filtered_in, and the
    corresponding phase value
    """
    # Min value
    if min_or_max == "min":
        g2_minmax = min(g2_initial_filtered_in)
    elif min_or_max == "max":
        g2_minmax = max(g2_initial_filtered_in)

    # Index where minimum is
    index = np.where(g2_initial_filtered_in == g2_minmax)[0][0]
    # print(index)
    
    # Phase value
    width_out = widths_in[index]
    
    # g2 value
    g2_out = g2_initial_filtered_in[index]
    
    print(("g2(0) = {:.5f} at width = {:.5f}"
           ).format(g2_out, width_out))
    
    return width_out, g2_out

#-----------------------------------------------------------------------------#
#                                  DATA THINGS                                #
#-----------------------------------------------------------------------------#
# Parameters
gamma, Omega, N, w0 = \
    np.genfromtxt(par_C_multi, delimiter="=", usecols=1, max_rows=7)
N = int(N)

#---------------------------------#
#     Read Data: Central Peak     #
#---------------------------------#
# Data: Multi-mode
K_C_multi  = np.genfromtxt(dat_C_multi, usecols=0, dtype='float')[1:]
g2_C_multi = np.genfromtxt(dat_C_multi, usecols=1, dtype='float')[1:]

# Data: Single-mode
K_C_single  = np.genfromtxt(dat_C_single, usecols=0, dtype='float')[1:]
g2_C_single = np.genfromtxt(dat_C_single, usecols=1, dtype='float')[1:]

#-------------------------------#
#     Read Data: Right Peak     #
#-------------------------------#
# Data: Multi-mode
K_R_multi  = np.genfromtxt(dat_R_multi, usecols=0, dtype='float')[1:]
g2_R_multi = np.genfromtxt(dat_R_multi, usecols=1, dtype='float')[1:]

# Data: Single-mode
K_R_single  = np.genfromtxt(dat_R_single, usecols=0, dtype='float')[1:]
g2_R_single = np.genfromtxt(dat_R_single, usecols=1, dtype='float')[1:]

#-----------------------------------------------------------------------------#
#                        PLOT g^{(2)} (SINGLE-FIGURE)                         #
#-----------------------------------------------------------------------------#
fig = plt.figure("g2 Initial Value")
ax = plt.gca()

# Plot filtered correlations
# Multi-mode
ax.semilogx(K_C_multi, g2_C_multi, color='C0', ls='solid', label='Multi-Mode') 
# Single-mode
ax.semilogx(K_C_single, g2_C_single, color='C1', ls='dashed', label='Single-Mode')
  
# Legend
ax.legend(loc='lower left')

# Labels
ax.set_xlabel(r'$K / \gamma$')
ax.set_ylabel(r'$g^{(2)}(0, 0; 0, 0)$')

# Ticks
ax.set_yticks(np.arange(0.0, 3.5, 0.5))
ax.set_yticks(np.arange(0.25, 3.5, 0.5), minor=True)

# Grid stuff
ax.grid()
    
# Axis lims
ax.set_ylim(-0.1, 2.70)
ax.set_xlim(8e-6, 1.5e2)

fig.tight_layout()
# fig.savefig(filename_out)
fig.show()

#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
fig, ax = plt.subplots(num='Combined Initial Scan', sharex=False, nrows=2, ncols=1,
                       figsize=[6, 6])

#----------------------#
#     Central Peak     #
#----------------------#
# Plot multi-mode
ax[0].semilogx(K_C_multi, g2_C_multi, color='C0', ls='solid',
               label=r'Multi-Mode')
# Plot single-mode
ax[0].semilogx(K_C_single, g2_C_single, color='C1', ls='dashed', alpha=0.85,
               label=r'Single-Mode')

#--------------------#
#     Right Peak     #
#--------------------#
# Plot multi-mode
ax[1].semilogx(K_R_multi, g2_R_multi, color='C0', ls='solid',
               label=r'Multi-Mode')
# Plot single-mode
ax[1].semilogx(K_R_single, g2_R_single, color='C1', ls='dashed', alpha=0.85,
               label=r'Single-Mode')

#------------------#
#     Add Text     #
#------------------#
# Central peak
ax[0].text(x=1.5e-6, y=2.55, s='(a)')
ax[1].text(x=1.5e-6, y=1.95, s='(b)')

#--------------------#
#     Axis Stuff     #
#--------------------#
# Cycle through axes
for j in range(2):
    # Legend
    ax[j].legend(loc='lower left')
    
    # Grid
    ax[j].grid()
    
    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    ax[j].set_yticks(np.arange(0.0, 3.5, 0.5))
    ax[j].set_yticks(np.arange(0.25, 3.5, 0.5), minor=True)

    #---------------------#
    #     Axis Limits     #
    #---------------------#
    # Axis lims
    ax[j].set_xlim(8e-6, 1.25e2)
    
    if j == 0:
        ax[j].set_ylim(-0.05, 2.70)
    elif j == 1:
        ax[j].set_ylim(-0.05, 2.1)

    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax[j].set_xlabel(r'$K / \gamma$')
    if j == 0:
        ax[j].set_ylabel(r'$g^{(2)}(0, 0; 0, 0)$')
    elif j == 1:
        ax[j].set_ylabel(r'$g^{(2)}(\Omega, 0; \Omega, 0)$')

#----------------------#
#     Figure Stuff     #
#----------------------#
# fig.tight_layout(pad=0)
fig.savefig(filename_out)
fig.show()