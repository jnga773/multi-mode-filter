#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 11:29:22 2024

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../../paper_style.mplstyle')

plt.close('all')

# Directories: Single Mode
dir_C_single = "./data_files/SingleMode/Centre/"
dir_R_single = "./data_files/SingleMode/Right/"
# Directories: Multi Mode
dir_C_multi = "./data_files/MultiMode/Centre/"
dir_R_multi = "./data_files/MultiMode/Right/"


# Figure filename
filename_out = "../../../images/sect5/fig10_g2_width_initial.pdf"

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def calc_min_initial(x_in, y_in, peak_in):
    """
    Calculates the minimum value of g2_initial_filtered_in, and the
    corresponding phase value
    """
    from numpy import where
    # Min value
    if peak_in == 'right':
        index = where(y_in == min(y_in))[0][0]
    elif peak_in == 'centre':
        index = where(abs(y_in - 1) == min(abs(y_in - 1)))[0][0]
        
    y_min = y_in[index]

    # Index where minimum is
    # print(index)
    
    # Phase value
    x_min = x_in[index]
    
    print(("Minimum g2(0) = {:.5f}, at width = {:.3f}"
           ).format(y_min, x_min))
    
    return x_min, y_min

#-----------------------------------------------------------------------------#
#                                  DATA THINGS                                #
#-----------------------------------------------------------------------------#
# Parameters
gamma, Omega, N, m, w0 = \
    np.genfromtxt(dir_C_multi + "parameters_initial.txt", usecols=1, delimiter='=')

# Read halfwidths
widths = np.genfromtxt(dir_C_single + "g2_comparison_average.txt",
                       skip_header=1, usecols=0)

#---------------------------------#
#     Read Data: Central-Peak     #
#---------------------------------#
# Multi-mode
g2_C_multi = np.genfromtxt(dir_C_multi + "g2_initial.txt",
                           skip_header=1, usecols=3)

# Single-Mode
g2_C_single = np.genfromtxt(dir_C_single + "g2_initial.txt",
                            skip_header=1, usecols=3)

#-------------------------------#
#     Read Data: Right-Peak     #
#-------------------------------#
# Multi-mode
g2_R_multi = np.genfromtxt(dir_R_multi + "g2_initial.txt",
                           skip_header=1, usecols=3)

# Single-Mode
g2_R_single = np.genfromtxt(dir_R_single + "g2_initial.txt",
                            skip_header=1, usecols=3)

#------------------------------#
#     Print Minimum Points     #
#------------------------------#
# Central Peak
print("Central Peak | Multi-Mode")
min_width_C_multi, min_g2_C_multi = \
    calc_min_initial(widths[0:100], g2_C_multi[0:100], 'centre')
print('')

print("Central Peak | Single-Mode")
min_width_C_single, min_g2_C_single = \
    calc_min_initial(widths[0:100], g2_C_single[0:100], 'centre')
print('')

# Right peak
print("Right Peak   | Multi-Mode")
min_width_R_multi, min_g2_R_multi = \
    calc_min_initial(widths[0:100], g2_R_multi[0:100], 'right')
print('')

print("Right Peak   | Single-Mode")
min_width_R_single, min_g2_R_single = \
    calc_min_initial(widths[0:100], g2_R_single[0:100], 'right')


#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
fig, ax = plt.subplots(num="g2 Difference (Width Scan)", nrows=2, ncols=1,
                       sharex=False, figsize=[6, 6])

#----------------------------#
#     Plot: Central Peak     #
#----------------------------#
# Plot: Single-Mode
ax[0].plot(widths, g2_C_multi, color='C0', ls='solid',
           label='Multi-Mode')
ax[0].plot(widths, g2_C_single, color='C1', ls='dashed',
           label='Single-Mode')

# Add crosses
# ax[0].plot(min_width_C_multi, min_g2_C_multi, marker='x', color='C0')
# ax[0].plot(min_width_C_single, min_g2_C_single, marker='x', color='C1')

# Add line for secular apprixmation initial value
ax[0].plot([0, max(widths)], [1, 1], color='k', ls='dotted')

# Legend
ax[0].legend(loc='upper left')

#--------------------------#
#     Plot: Right-Peak     #
#--------------------------#
ax[1].plot(widths, g2_R_multi, color='C0', ls='solid',
           label='Multi-Mode')
ax[1].plot(widths, g2_R_single, color='C1', ls='dashed',
           label='Single-Mode')

# Add crosses
# ax[1].plot(min_width_R_multi, min_g2_R_multi, marker='x', color='C0')
# ax[1].plot(min_width_R_single, min_g2_R_single, marker='x', color='C1')

# Add line for secular apprixmation initial value
ax[1].plot([0, max(widths)], [0, 0], color='k', ls='dotted')

# Legend
ax[1].legend(loc='upper left')

#-------------------------#
#     Add Figure Text     #
#-------------------------#
ax[0].text(x=-2.75, y=2.00, s='(a)')
ax[1].text(x=-2.75, y=1.75, s='(b)')

#--------------------#
#     Axis Ticks     #
#--------------------#
# Ticks: X-axis
ax[0].set_xticks(np.arange(0.0, 25.0, 2.5))
ax[0].set_xticks(np.arange(1.25, 25, 2.5), minor=True)

ax[1].set_xticks(np.arange(0.0, 25.0, 2.5))
ax[1].set_xticks(np.arange(1.25, 25, 2.5), minor=True)

# Ticks: Y-axis
ax[0].set_yticks(np.arange(0.0, 2.5, 0.25))
ax[0].set_yticks(np.arange(0.125, 2.5, 0.25), minor=True)

ax[1].set_yticks(np.arange(0.0, 2.0, 0.25))
ax[1].set_yticks(np.arange(0.125, 2.0, 0.25), minor=True)

# Tick parameters
ax[0].tick_params(which='minor', axis='both', grid_linestyle='dashed')
ax[1].tick_params(which='minor', axis='both', grid_linestyle='dashed')

#---------------------#
#     Axis Limits     #
#---------------------#
ax[0].set_xlim(-0.2, 20.2)
ax[0].set_ylim(-0.05, 2.05)

ax[1].set_xlim(-0.2, 20.2)
ax[1].set_ylim(-0.05, 1.8)

#---------------#
#     Grids     #
#---------------#
ax[0].grid(which='both')
ax[1].grid(which='both')

#---------------------#
#     Axis Labels     #
#---------------------#
ax[0].set_xlabel(r'$K / \gamma$')
ax[0].set_ylabel(r'$g^{(2)}(0, 0; 0, 0)$')

ax[1].set_xlabel(r'$K / \gamma$')
ax[1].set_ylabel(r'$g^{(2)}(\Omega, 0; \Omega, 0)$')

#----------------------#
#     Figure Stuff     #
#----------------------#
fig.tight_layout()
fig.savefig(filename_out)
fig.show()
