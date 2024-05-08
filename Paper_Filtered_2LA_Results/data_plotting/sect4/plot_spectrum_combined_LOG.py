#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 13:00:59 2023

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../paper_style.mplstyle')

plt.close('all')

# Figure filename
filename_out = "../../images/sect4/fig8_combined_g1_filtered.pdf"

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def spectrum(tau_input, corr_input):
    from numpy.fft import fft, fftshift, fftfreq
    from numpy import where, mean, pi, array
    from numpy import max as npmax
    
    # Shift the arrays so they are arranged from negative to positive freq
    fft = fft(corr_input)#, norm='ortho')
    fft = fftshift(fft)
    freq = fftfreq(tau_input.shape[0], tau_input[1]-tau_input[0])
    freq = fftshift(freq)
    
    # As the central peak is a delta function, we ignore the plot for w=0. To
    # do this we find the index in wlist where it is equal to 0.0, ignore it,
    # and create a new list excluding this value.
    indices = where(freq != 0.0)[0]
        
    # Remove zero frequency term
    spec_output = fft[indices]
    
    # Take only the real part
    spec_output = spec_output.real
    
    # take away non zero tails
    spec_output = spec_output - mean(spec_output[0])
    wlist_output = freq[indices] # wlist is in terms of FFT frequencies
    wlist_output = array(wlist_output) * 2 * pi
    
    # Normalise
    spec_output = spec_output / npmax(spec_output)
    
    return spec_output, wlist_output

def calc_mollow_triplet(gamma_in, Omega_in, tau_in):
    """
    Calculates the power spectrum (Mollow triplet) for a driven two-level atom.
    """
    from qutip import sigmam, correlation_2op_1t
    from numpy import sqrt
    
    # Operators
    sm = sigmam()
    
    # Hamiltonian
    H = 0.5 * Omega_in * (sm + sm.dag())
    
    # Collapse and evaluation operators
    c_ops = [sqrt(gamma_in) * sm]
    
    #----------------------------------------------------#
    #     Calculate First-Order Correlation Function     #
    #----------------------------------------------------#
    # Calculate the first-order correlation function
    G1_out = correlation_2op_1t(H, None, tau_in, c_ops, a_op=sm.dag(), b_op=sm,
                                reverse=False)
    # Calculate spectrum
    mollow_out, wlist_f = spectrum(tau_in, G1_out)
    
    return mollow_out.real, wlist_f

#-----------------------------------------------------------------------------#
#                                  DATA THINGS                                #
#-----------------------------------------------------------------------------#
# Read parameters
gamma, Omega, N, kappa, dw, w0, dt, tau_max = \
    np.genfromtxt("./data_files/right/g1_parameters_multi.txt", delimiter="=", usecols=1,
                  max_rows=8)
N = int(N)

# Calculate Mollow triplet
spec_mollow, wlist_mollow = calc_mollow_triplet(gamma, Omega, np.arange(0, tau_max + dt, dt))

# Pull halfwidth values
halfwidths = [1, 2, 4, 8, 16]

halfwidth_plot_indices = [1, 3, 4]

# Read frequencies
wlist = np.genfromtxt('./data_files/right/spec_single.txt', usecols=0, dtype='float')

# Read data
spec_single = []
spec_multi = []

for peak in ['centre', 'right']:
    # Single-mode data
    spec = np.genfromtxt('./data_files/{}/spec_single.txt'.format(peak), dtype='float',
                         usecols=(1, 2, 3, 4, 5))
    # Append to array
    spec_single.append(spec)

    # Multi-mode data
    spec = np.genfromtxt('./data_files/{}/spec_multi.txt'.format(peak), dtype='float',
                         usecols=(1, 2, 3, 4, 5))
    spec_multi.append(spec)

#-----------------------------------------------------------------------------#
#                               PLOT SPECTRUM                                 #
#-----------------------------------------------------------------------------#
# Plot colours
plot_colours = ['C0', 'C1', 'C2']
# Plot linestyles
plot_linestyle = ['solid', 'dashed', 'dashdot']

# Create figure
fig, ax = plt.subplots(num='Spectra combined', nrows=2, ncols=2, figsize=[12, 7],
                       sharex=False, sharey=False)
        
#---------------------------#
#     Plot: Single-Mode     #
#---------------------------#
# for i in range(len(halfwidths)):
for i in range(len(halfwidth_plot_indices)):
    # Data index
    i_data = halfwidth_plot_indices[i]
    # Plot colour
    colour = plot_colours[i]
    # Plot linestyle
    linestyle = plot_linestyle[i]
    
    # CENTRAL PEAK
    # Grab data:
    spec = spec_single[0]
    spec_plot = spec[:, i_data]
    halfwidth = halfwidths[i_data]
    
    # Plot
    # label_str = rf'Single-Mode $\left( K / \gamma = {halfwidth} \right)$'
    label_str = rf'Single-Mode $\left( K = {halfwidth} \gamma \right)$'
    ax[0, 0].semilogy(wlist, spec_plot, color=colour, ls=linestyle,
                      label=label_str)
    
    # Right PEAK
    # Grab data:
    spec = spec_single[1]
    spec_plot = spec[:, i_data]
    halfwidth = halfwidths[i_data]
    
    # Plot
    # label_str = rf'Single-Mode $\left( K / \gamma = {halfwidth} \right)$'
    label_str = rf'Single-Mode $\left( K = {halfwidth} \gamma \right)$'
    ax[0, 1].semilogy(wlist, spec_plot, color=colour, ls=linestyle,
                      label=label_str)

#--------------------------#
#     Plot: Multi-Mode     #
#--------------------------#
# for i in range(len(halfwidths)):
for i in range(len(halfwidth_plot_indices)):
    # Data index
    i_data = halfwidth_plot_indices[i]
    # Plot colour
    colour = plot_colours[i]
    # Plot linestyle
    linestyle = plot_linestyle[i]
    
    # CENTRAL PEAK
    # Grab data:
    spec = spec_multi[0]
    spec_plot = spec[:, i_data]
    halfwidth = halfwidths[i_data]
    
    # Plot
    # label_str = rf'Multi-Mode $\left( K / \gamma = {halfwidth} \right)$'
    label_str = rf'Multi-Mode $\left( K = {halfwidth} \gamma \right)$'
    ax[1, 0].semilogy(wlist, spec_plot, color=colour, ls=linestyle,
                      label=label_str)
    
    # Right PEAK
    # Grab data:
    spec = spec_multi[1]
    spec_plot = spec[:, i_data]
    halfwidth = halfwidths[i_data]
    
    # Plot
    # label_str = rf'Multi-Mode $\left( K / \gamma = {halfwidth} \right)$'
    label_str = rf'Multi-Mode $\left( K = {halfwidth} \gamma \right)$'
    ax[1, 1].semilogy(wlist, spec_plot, color=colour, ls=linestyle,
                      label=label_str)
    
#--------------------------#
#     Plot: Unfiltered     #
#--------------------------#
for i in range(2):
    for j in range(2):
        # Plot Mollow triplet
        ax[i, j].plot(wlist_mollow, spec_mollow, color='k', alpha=0.75,
                      ls='dotted', label='Mollow Triplet')

#------------------#
#     Add Text     #
#------------------#
# Single-mode
ax[0, 0].text(x=-5.5, y=4.8, s='(a)')
ax[0, 1].text(x=-5.5, y=4.8, s='(b)')
# Multi-mode
ax[1, 0].text(x=-5.5, y=4.8, s='(c)')
ax[1, 1].text(x=-5.5, y=4.8, s='(d)')

#--------------------#
#     Axis Stuff     #
#--------------------#
# Cycle through axes
for i in range(2):
    for j in range(2):
        # Legend
        ax[i, j].legend()
        
        # Grid
        ax[i, j].grid()
        
        #--------------------#
        #     Axis Ticks     #
        #--------------------#
        ax[i, j].set_xticks(np.arange(-2*np.pi, 7*np.pi, np.pi))
        ax[i, j].set_xticks(np.arange(-1.5*np.pi, 6*np.pi, np.pi), minor=True)
        ax[i, j].set_xticklabels([r'$-2\pi$', r'$-\pi$', '0', r'$\pi$', r'$2 \pi$', r'$3 \pi$', r'$4 \pi$', r'$5 \pi$', r'$6 \pi$'])

        # ax[i, j].set_yticks(np.arange(0.0, 1.2, 0.2))
        # ax[i, j].set_yticks(np.arange(0.1, 1.1, 0.2), minor=True)

        #---------------------#
        #     Axis Limits     #
        #---------------------#
        ax[i, j].set_xlim(-np.pi, 6 * np.pi)
        ax[i, j].set_ylim(1e-4, 1e1)

        #---------------------#
        #     Axis Labels     #
        #---------------------#
        ax[i, j].set_xlabel(r'$( \omega - \omega_{A} ) / \gamma$')
        # ax.set_ylabel(r'Power Spectrum (a.u.)')
        ax[i, j].set_ylabel(r'$S_{\mathrm{inc}}(\omega)$ (a.u.)')

#----------------------#
#     Figure Stuff     #
#----------------------#
# fig.tight_layout(pad=0)
fig.savefig(filename_out)
fig.show()
