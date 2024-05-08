#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 13:00:05 2023

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../paper_style.mplstyle')

plt.close('all')

# Filename for figure
filename_out = "../../images/sect2/fig4_multi_mode_frequency_response.pdf"

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def calc_alpha(w_in, w0_in, kappa_in, E_in):
    """
    Calculates the long-time limit coherent field amplitude:
    \alpha(t) = E * e^{-i \omega t} / (\kappa + i(\omega_{0} - \omega)t)
    """
    # Big time for long time limit (lol)
    t_f = 1E6

    # Semi-classical field amplitude
    alpha_f = E_in / (kappa_in + 1j * (w0_in - w_in))
    alpha_f *= np.exp(-1j * w_in * t_f)

    return alpha_f

def multi_mode_alpha(omega_in, w0_in, kappa_in, E0_in, N_in, dw_in, m_in):
    """
    Calculates the collective mode alpha by summing over all modes
    """
    # Array of mode numbers
    modes = np.arange(-N_in, N_in + 1, 1)

    # Create empty list for A_out to be added to
    A_out = np.zeros(len(omega_in), dtype='complex128')

    # cycle through modes and calculate wj and epsj
    for j in modes:
        # Resonance frequency
        wj = w0_in + (j * dw_in)

        # Coupling constant
        Ej = (E0_in / np.sqrt(2*N_in + 1)) * np.exp(1j * m_in * j * np.pi / N_in)
        # Ej = eps0_in * np.cos(j * np.pi / N_in)

        # Decay rate
        kappaj = kappa_in
        # kappaj = kappa_in * np.cos(j + 0.5 * np.pi / N_in)

        # Calculate and add the cavity amplitude
        A_out += calc_alpha(omega_in, wj, kappaj, Ej)

    return A_out

def filter_frequency_response(omega_in, w0_in, kappa_in, E0_in, N_in, dw_in, m_in):
    """
    Calculate the frequency response of the filter by taking the mode square
    """
    # Calculate the cavity amplitude
    if N_in == 0:
        alpha = calc_alpha(omega_in, w0_in, kappa_in, E0_in)
    else:
        alpha = multi_mode_alpha(omega_in, w0_in, kappa_in, E0_in, N_in,
                                 dw_in, m_in)

    # Calculate the response
    response = np.abs(alpha) ** 2
    # Normalise to peak value
    response *= 1 / max(response)

    return response

#-----------------------------------------------------------------------------#
#                                 PARAMETERS                                  #
#-----------------------------------------------------------------------------#
# Number of modes either side of central mode (2N+1 total modes)
# N = 40
N = 80

# Mode spacing
dw1 = 1.0 / N
dw2 = 2.0 / N
dw3 = 4.0 / N

# Cavity decay rate - multi-mode
# kappa = 0.2
kappa = 0.07

# Central mode frequency
w0 = 0.0

# Driving amplitude
E0 = 1.0

# Phase multiplier
m0 = 0
m1 = 1
m2 = 2

# Array of frequencies to calculate over
omega = np.arange(w0-10, w0+10.0, 0.001)

#-----------------------------------------------------------------------------#
#                             CALCULATE REPSONSES                             #
#-----------------------------------------------------------------------------#

# Single mode response
# a_single = alpha(omega, w0, kappa_single, eps0)
single_response = filter_frequency_response(omega, w0, 1.0, E0, 0, 0, 0)

# Multi-mode response: Different bandwidths
response_dw1 = [filter_frequency_response(omega, w0, kappa, E0, N, dw1, m0),
                filter_frequency_response(omega, w0, kappa, E0, N, dw1, m1),
                filter_frequency_response(omega, w0, kappa, E0, N, dw1, m2)]

response_dw2 = [filter_frequency_response(omega, w0, kappa, E0, N, dw2, m0),
                filter_frequency_response(omega, w0, kappa, E0, N, dw2, m1),
                filter_frequency_response(omega, w0, kappa, E0, N, dw2, m2)]

response_dw3 = [filter_frequency_response(omega, w0, kappa, E0, N, dw3, m0),
                filter_frequency_response(omega, w0, kappa, E0, N, dw3, m1),
                filter_frequency_response(omega, w0, kappa, E0, N, dw3, m2)]

#-----------------------------------------------------------------------------#
#                                    PLOT                                     #
#-----------------------------------------------------------------------------#
# fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=False,
#                        figsize=[6, 7.5])
fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=False,
                       figsize=[6, 7])


# Single-mode
for i in range(3):
    #---------------------------#
    #     Plot: Single-Mode     #
    #---------------------------#
    # label_str = r'$N = {}, K / \left| \mathcal{{E}}_{{d}} \right| = {}$'.format(0, 1)
    label_str = r'$N = {}, K = \left| \mathcal{{E}}_{{d}} \right|$'.format(0)
    ax[i].plot(omega, single_response, color='k', ls='dotted', alpha=0.75,
               label=label_str)

    #--------------------------#
    #     Plot: Multi-Mode     #
    #--------------------------#
    # Plot N \delta\omega1
    # label_str = r'$N \delta\omega / \left| \mathcal{{E}}_{{d}} \right| = {}$'.format(int(N * dw1)))
    # label_str = r'$N = {}, K / \left| \mathcal{{E}}_{{d}} \right| = {}$'.format(N, int(N * dw1))
    label_str = r'$N = {}, K = \left| \mathcal{{E}}_{{d}} \right|$'.format(N)
    ax[i].plot(omega, response_dw1[i], color='C0', ls='solid', 
              label=label_str)
    
    # Plot N \delta\omega2
    # label_str = r'$N \delta\omega / \left| \mathcal{{E}}_{{d}} \right| = {}$'.format(int(N * dw2)))
    # label_str = r'$N = {}, K / \left| \mathcal{{E}}_{{d}} \right| = {}$'.format(N, int(N * dw2))
    label_str = r'$N = {}, K = {} \left| \mathcal{{E}}_{{d}} \right|$'.format(N, int(N * dw2))
    ax[i].plot(omega, response_dw2[i], color='C1', ls='dashed',
               label=label_str)
    
    # Plot N \delta\omega3
    # label_str = r'$N \delta\omega / \left| \mathcal{{E}}_{{d}} \right| = {}$'.format(int(N * dw3)))
    # label_str = r'$N = {}, K / \left| \mathcal{{E}}_{{d}} \right| = {}$'.format(N, int(N * dw3))
    label_str = r'$N = {}, K = {} \left| \mathcal{{E}}_{{d}} \right|$'.format(N, int(N * dw3))
    ax[i].plot(omega, response_dw3[i], color='C2', ls='dashdot',
               label=label_str)
    
    #----------------#
    #     Legend     #
    #----------------#
    ax[i].legend()

    #--------------------#
    #     Axis Ticks     #
    #--------------------#
    # ax[i].set_xticks(np.arange(w0-10, w0+12, 2), minor=False)
    # ax[i].set_xticks(np.arange(w0-9, w0+10, 2), minor=True)
    ax[i].set_xticks(np.arange(-5, 6, 1))
    ax[i].set_xticks(np.arange(-4.5, 6, 1), minor=True)
    
    ax[i].set_yticks(np.arange(0.0, 1.2, 0.2))
    ax[i].set_yticks(np.arange(0.1, 1.3, 0.2), minor=True)

    #---------------------#
    #     Axis limits     #
    #---------------------#
    # ax.set_xlim(-10.5, 10.5)
    ax[i].set_xlim(-5, 5)
    ax[i].set_ylim(-0.025, 1.025)
    # ax[i].set_ylim(-0.05, 1.1)
    
    #---------------------#
    #     Axis Labels     #
    #---------------------#
    # ax.set_xlabel(r'$\left( \omega - \omega_{c} \right) / |\mathcal{E}_{d}|$')
    ax[i].set_ylabel(r'$\left| \bar{A}(\omega) \right|^{2}$ (a.u.)')
    
    #----------------------#
    #     Figure stuff     #
    #----------------------#
    # Grid
    ax[i].grid()

# Add text for sub labels
ax[0].text(x=-4.8, y=0.9, s='(a)')
ax[1].text(x=-4.8, y=0.9, s='(b)')
ax[2].text(x=-4.8, y=0.9, s='(c)')
    
# y-Axis label
ax[2].set_xlabel(r'$\omega / |\mathcal{E}_{d}|$')

#----------------------------------------#
#     Arrows Labelling Lines (m = 0)     #
#----------------------------------------#
# Arrow properties
arrow_props = {'arrowstyle': '->', 'facecolor': 'k'}

#----------------------------------------#
#     Arrows Labelling Lines (m = 1)     #
#----------------------------------------#

#----------------------------------------#
#     Arrows Labelling Lines (m = 2)     #
#----------------------------------------#


# fig.tight_layout(pad=0)
fig.savefig(filename_out)
fig.show()
