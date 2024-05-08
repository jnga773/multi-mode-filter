# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 12:40:36 2024

@author: nexus
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../../paper_style.mplstyle')

plt.close('all')

# Filenames: Right peak
dat_R_multi  = "./data_files/right/g2_corr_multi.txt"
dat_R_single = "./data_files/right/g2_corr_single.txt"

par_R_multi  = "./data_files/right/g2_parameters_multi.txt"
par_R_single = "./data_files/right/g2_parameters_single.txt"

# Filenames: Central peak
dat_C_multi  = "./data_files/centre/g2_corr_multi.txt"
dat_C_single = "./data_files/centre/g2_corr_single.txt"

par_C_multi  = "./data_files/centre/g2_parameters_multi.txt"
par_C_single = "./data_files/centre/g2_parameters_single.txt"

# Figure filename
filename_out = '../../../images/sect5/fig9_combined_g2_width_scan.pdf'

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
gamma, Omega, N, kappa, w0, dt, tau_max = \
    np.genfromtxt(par_R_multi, delimiter="=", usecols=1, max_rows=7)
N = int(N)

# Halfwidths
halfwidths = np.genfromtxt(par_R_multi, skip_header=9, usecols=0)
dw         = np.genfromtxt(par_R_multi, skip_header=9, usecols=1)
kappa      = np.genfromtxt(par_R_multi, skip_header=9, usecols=2)

# Read tau
tau = np.arange(0.0, tau_max + dt, dt)

# Dressed state correlation
g2_R = g2_dressed_states(tau, Omega, gamma, Omega)
g2_C = g2_dressed_states(tau, Omega, gamma, 0.0)

# Single-mode data
g2_R_single = np.genfromtxt(dat_R_single)
g2_C_single = np.genfromtxt(dat_C_single)

# Multi-mode data
g2_R_multi = np.genfromtxt(dat_R_multi)
g2_C_multi = np.genfromtxt(dat_C_multi)

# Halfwidths that I will plot
halfwidths_plot = [2, 4, 8]

#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
# Plot colours
plot_colours = ['C0', 'C1', 'C2']
# Plot linestyles
plot_linestyle = ['solid', 'dashed', 'dashdot']


fig, ax = plt.subplots(num='Good vs Bad', sharex=False, nrows=2, ncols=2,
                       figsize=[12, 7])

#----------------------------------#
#     Plot: Cycle through Data     #
#----------------------------------#
# for i in range(len(halfwidths)):
for i in range(len(halfwidths_plot)):
    # Halfwidth for plot
    K_plot = halfwidths_plot[i]
    # Plot colour
    colour = plot_colours[i]
    # Plot linestyle
    linestyle = plot_linestyle[i]
    
    #---------------------------#
    #     Plot: Single-Mode     #
    #---------------------------#
    # CENTRAL PEAK
    # Grab data:
    g2_plot = g2_C_single[:, K_plot-1]
    
    # Plot
    # label_str = rf'Single-Mode $\left( K / \gamma = {K_plot} \right)$'
    label_str = rf'Single-Mode $\left( K = {K_plot} \gamma \right)$'
    ax[0, 0].plot(tau, g2_plot, color=colour, ls=linestyle,
                  label=label_str)
    
    # RIGHT PEAK
    # Grab data:
    g2_plot = g2_R_single[:, K_plot-1]
    
    # Plot
    # label_str = rf'Single-Mode $\left( K / \gamma = {K_plot} \right)$'
    label_str = rf'Single-Mode $\left( K = {K_plot} \gamma \right)$'
    ax[0, 1].plot(tau, g2_plot, color=colour, ls=linestyle,
                  label=label_str)
    
    #--------------------------#
    #     Plot: Multi-Mode     #
    #--------------------------#    
    # CENTRAL PEAK
    # Grab data:
    g2_plot = g2_C_multi[:, K_plot-1]
    
    # Plot
    # label_str = rf'Multi-Mode $\left( K / \gamma = {K_plot} \right)$'
    label_str = rf'Multi-Mode $\left( K = {K_plot} \gamma \right)$'
    ax[1, 0].plot(tau, g2_plot, color=colour, ls=linestyle,
                  label=label_str)
    
    # RIGHT PEAK
    # Grab data:
    g2_plot = g2_R_multi[:, K_plot-1]
    
    # Plot
    # label_str = rf'Multi-Mode $\left( K / \gamma = {K_plot} \right)$'
    label_str = rf'Multi-Mode $\left( K = {K_plot} \gamma \right)$'
    ax[1, 1].plot(tau, g2_plot, color=colour, ls=linestyle,
                  label=label_str)
    
#--------------------------#
#     Plot: Unfiltered     #
#--------------------------#
for i in range(2):
    for j in range(2):
        if j == 0:
            # Plot central peak
            ax[i, j].plot(tau, g2_C, color='k', alpha=0.75,
                          ls='dotted', label='Secular Approximation Eq. (38a)')
        elif j == 1:
            # Plot right peak
            ax[i, j].plot(tau, g2_R, color='k', alpha=0.75,
                          ls='dotted', label='Secular Approximation Eq. (38c)')


#------------------#
#     Add Text     #
#------------------#
# Central peak
ax[0, 0].text(x=-0.65, y=1.3, s='(a)')
ax[0, 1].text(x=-0.65, y=1.0, s='(b)')

# Right peak
ax[1, 0].text(x=-0.65, y=1.3, s='(c)')
ax[1, 1].text(x=-0.65, y=1.0, s='(d)')

#--------------------#
#     Axis Stuff     #
#--------------------#
# Cycle through axes
for i in range(2):
    for j in range(2):
        # Legend
        ax[i, j].legend(loc='lower right')
        
        # Grid
        ax[i, j].grid()
        
        #--------------------#
        #     Axis Ticks     #
        #--------------------#
        # ax[i, j].set_xticks(np.arange(0.0, 12, 2))
        # ax[i, j].set_xticks(np.arange(1.0, 11, 2), minor=True)
        ax[i, j].set_xticks(np.arange(0.0, 5.5, 1))
        ax[i, j].set_xticks(np.arange(0.5, 5.5, 1), minor=True)

        if j == 0:
            ax[i, j].set_yticks(np.arange(0.7, 1.5, 0.15))
            ax[i, j].set_yticks(np.arange(0.775, 1.5, 0.15), minor=True)
        elif j == 1:
            ax[i, j].set_yticks(np.arange(0.0, 1.25, 0.25))
            ax[i, j].set_yticks(np.arange(0.125, 1.25, 0.25), minor=True)

        # ax[i, j].set_yticks(np.arange(0.0, 1.2, 0.2))
        # ax[i, j].set_yticks(np.arange(0.1, 1.1, 0.2), minor=True)

        #---------------------#
        #     Axis Limits     #
        #---------------------#
        # ax[i, j].set_xlim(-0.1, 10.1)
        ax[i, j].set_xlim(-0.05, 5.05)

        if j == 0:
            ax[i, j].set_ylim(0.67, 1.33)
        elif j == 1:
            ax[i, j].set_ylim(-0.05, 1.05)

        #---------------------#
        #     Axis Labels     #
        #---------------------#
        ax[i, j].set_xlabel(r'$\gamma \tau$')
        if j == 0:
            ax[i, j].set_ylabel(r'$g^{(2)}(0, 0; 0, \tau)$')
        elif j == 1:
            ax[i, j].set_ylabel(r'$g^{(2)}(\Omega, 0; \Omega, \tau)$')

#----------------------#
#     Figure Stuff     #
#----------------------#
# fig.tight_layout(pad=0)
# fig.savefig(filename_out)
fig.show()