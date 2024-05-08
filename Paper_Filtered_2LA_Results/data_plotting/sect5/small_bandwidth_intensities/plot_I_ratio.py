#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:22:12 2024

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
par_C_multi  = './data_files/centre/intensity_ratio_parameters_multi.txt'
# # Parameters: Single-mode
# par_C_single = './data_files/centre/intensity_ratio_parameters_single.txt'

# Data: Multi-mode
dat_C_multi  = './data_files/centre/intensity_ratio_multi.txt'
# Data: Single-mode
dat_C_single = './data_files/centre/intensity_ratio_single.txt'

#--------------------#
#     Right Peak     #
#--------------------#
# # # Parameters: Multi-mode
# # par_R_multi  = './data_files/right/intensity_ratio_parameters_multi.txt'
# # # Parameters: Single-mode
# # par_R_single = './data_files/right/intensity_ratio_parameters_single.txt'

# # Data: Multi-mode
# dat_R_multi  = './data_files/right/intensity_ratio_multi.txt'
# # Data: Single-mode
# dat_R_single = './data_files/right/intensity_ratio_single.txt'

#-------------------------#
#     Figure Filename     #
#-------------------------#
filename_out = '../../../images/sect5/fig13_intensity_ratio.pdf'

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
K_C_multi = np.genfromtxt(dat_C_multi, usecols=0, dtype='float')[1:]
I_C_multi = np.genfromtxt(dat_C_multi, usecols=1, dtype='float')[1:]

# Data: Single-mode
K_C_single = np.genfromtxt(dat_C_single, usecols=0, dtype='float')[1:]
I_C_single = np.genfromtxt(dat_C_single, usecols=1, dtype='float')[1:]

#-------------------------------#
#     Read Data: Right Peak     #
#-------------------------------#
# # Data: Multi-mode
# K_R_multi = np.genfromtxt(dat_R_multi, usecols=0, dtype='float')[1:]
# I_R_multi = np.genfromtxt(dat_R_multi, usecols=1, dtype='float')[1:]

# # Data: Single-mode
# K_R_single = np.genfromtxt(dat_R_single, usecols=0, dtype='float')[1:]
# I_R_single = np.genfromtxt(dat_R_single, usecols=1, dtype='float')[1:]
    
#-----------------------------------------------------------------------------#
#                                PLOT RATIO                                   #
#-----------------------------------------------------------------------------#
fig = plt.figure(num='Combined Initial Scan')
ax = plt.gca()

#----------------------#
#     Central Peak     #
#----------------------#
# Plot multi-mode
ax.loglog(K_C_multi, I_C_multi, color='C0', ls='solid',
          label=r'Multi-Mode')
# Plot single-mode
ax.loglog(K_C_single, I_C_single, color='C1', ls='dashed', alpha=0.85,
          label=r'Single-Mode')

# Legend
ax.legend(loc='upper left')

# Grid
ax.grid()

#---------------#
#     Inset     #
#---------------#
# Set up inst axis
axins = ax.inset_axes([0.65, 0.15, 0.33, 0.5])

axins.plot(K_C_multi, I_C_multi, ls='solid', color='C0',
           label='Multi-Mode')
axins.plot(K_C_single, I_C_single, ls='dashed', color='C1',
           label='Single-Mode')

# Axis ticks
axins.set_xticks(np.arange(5.0, 25.0, 5.0), minor=False)
axins.set_xticks(np.arange(7.5, 25.0, 5.0), minor=True)

axins.set_yticks(np.arange(150.0, 550.0, 100.0), minor=False)
axins.set_yticks(np.arange(200.0, 550.0, 100.0), minor=True)

# Axis limits
axins.set_xlim(5.0, 20.0)
axins.set_ylim(150, 450)

ax.indicate_inset_zoom(axins, edgecolor="black")
    
#--------------------#
#     Axis Ticks     #
#--------------------#
# ax[j].set_yticks(np.arange(0.0, 3.5, 0.5))
# ax[j].set_yticks(np.arange(0.25, 3.5, 0.5), minor=True)

#---------------------#
#     Axis Limits     #
#---------------------#
# Axis lims
# ax.set_xlim(1e-5, 1e2)
# ax.set_ylim(1e-3, 1e3)

ax.set_xlim(8e-6, 1.5e2)
ax.set_ylim(7e-4, 1.5e3)

#---------------------#
#     Axis Labels     #
#---------------------#
ax.set_xlabel(r'$K / \gamma$')
ax.set_ylabel(r'$I_{\mathrm{inc}} / I_{\mathrm{coh}}$')

#----------------------#
#     Figure Stuff     #
#----------------------#
# fig.tight_layout(pad=0)
fig.savefig(filename_out)
fig.show()