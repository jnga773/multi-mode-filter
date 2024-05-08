#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 12:29:02 2020

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../paper_style.mplstyle')

plt.close('all')

# Filename for figure
filename_out = "../../images/sect2/fig3_impulse_driving_cavity_amplitude.pdf"

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def theta(t_in):
    """
    Heaviside function
    """
    from numpy import heaviside
    return heaviside(t_in, 0.0)

def calc_alpha(t_in, wj_in, kappa_in, Ej_in, phi_in):
    """
    Calculates the long-time limit coherent field amplitude:
    \alpha_{j}(t) = E_{j} * e^{-(kappa + 1j * wj) * t} * theta(t)
    """
    # Cavity amplitude in the long time limit
    alpha_out = Ej_in
    alpha_out *= np.exp(-(kappa_in + (1j * wj_in)) * t_in)
    alpha_out *= theta(t_in)
    
    return alpha_out

def multi_mode_alpha(t_in, w0_in, kappa_in, E0_in, N_in, dw_in, m_in, phi_in):
    """
    Calculates the collective mode alpha by summing over all modes
    """
    # Array of mode numbers
    modes = np.arange(-N_in, N_in + 1, 1)
    
    # Create empty list for A_out to be added to
    A_out = np.zeros(len(t_in), dtype='complex128')
    
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
        A_out += calc_alpha(t_in, wj, kappaj, Ej, phi_in)

    return A_out

def multi_analytic(t_in, w0_in, kappa_in, E0_in, N_in, dw_in, m_in, phi_in):
    # Argument of the sinc function
    x = (N_in * dw_in * t_in) - (m_in * np.pi)
    x *= 1 / np.pi
    
    # Calculate response
    A_out = 0j
    # A_out += 2 * N_in * E0_in / np.sqrt(2 * N_in + 1)
    A_out += 2 * N_in * E0_in / np.sqrt(2 * N_in + 1)
    A_out *= theta(t_in - phi_in)
    A_out *= np.exp(-(kappa_in + 1j * w0_in) * t_in)
    A_out *= np.sinc(x)
    
    # Normalise
    A_out *= 1 / np.sqrt((2 * N_in) + 1)
    
    return A_out.real

def sinc_response(t_in, N_in, dw_in, m_in, phi_in):
    # Calculates just a sinc function to compare with
    x = (N_in * dw_in * t_in) - (m_in * np.pi)
    x *= 1 / np.pi
    
    sinc_out = np.sinc(x) * theta(t_in - phi_in)
    
    return sinc_out.real

#-----------------------------------------------------------------------------#
#                                 PARAMETERS                                  #
#-----------------------------------------------------------------------------#
# Number of modes either side of central mode (2N+1 total modes)
N = 20

# Mode spacing
# dw = 1.0 / N
dw = 2.0 / N
# dw = 4.0 / N

# Cavity decay rate - multi-mode
kappa = 0.2

# Central mode frequency
w0 = 0.0

# Driving amplitude
E0 = 1.0

# Phase multiplier
m0 = 0
m1 = 1
m2 = 2

# Time delay
phi = 0.0
# phi = 10.0

# Array of frequencies to calculate over
dt = 0.001
t_max = 3 * np.pi
t = np.arange(0.0, t_max + dt, dt)

#-----------------------------------------------------------------------------#
#                             CALCULATE REPSONSES                             #
#-----------------------------------------------------------------------------#
# Perfect sinc response
shifted_sinc = sinc_response(t, N, dw, m1, phi)
# shifted_sinc = multi_analytic(t, 0.0, 0.0, E0, N, dw, m1, phi)

# Multi mode response (with constant phase coupling m = 0)
# response_m0 = multi_mode_alpha(t, w0, kappa, E0, N, dw, m0, phi)
response_m0 = multi_analytic(t, w0, kappa, E0, N, dw, m0, phi)

# Multi mode response (with mode-dependent phase coupling m = 1)
# response_m1 = multi_mode_alpha(t, w0, kappa, E0, N, dw, m1, phi)
response_m1 = multi_analytic(t, w0, kappa, E0, N, dw, m1, phi)

# Multi mode response (with mode-dependent phase coupling m = 2)
# response_m2 = multi_mode_alpha(t, w0, kappa, E0, N, dw, m2, phi)
response_m2 = multi_analytic(t, w0, kappa, E0, N, dw, m2, phi)

#--------------------------------#
#     Normalisation constant     #
#--------------------------------#
# # Normalisation constant is 2N / sqrt(2N + 1)
# norm = (2 * N) / np.sqrt((2 * N) + 1)

# # Normalise all responses by this value
# response_m0  *= 1 / norm
# response_m1  *= 1 / norm
# response_m2  *= 1 / norm

#-----------------------------------------------------------------------------#
#                              PLOT: TIME-SERIES                              #
#-----------------------------------------------------------------------------#
fig = plt.figure(num='Impulse Driving')
ax = plt.gca()

#--------------#
#     Plot     #
#--------------#
# Multi_mode: m = 0
ax.plot(t, response_m0, color='C0', ls='solid',
        label=r'$m = {}$'.format(m0))

# Multi_mode: m = 1
ax.plot(t, response_m1, color='C1', ls='dashed',
        label=r'$m = {}$'.format(m1))

# Multi_mode: m = 2
ax.plot(t, response_m2, color='C2', ls='dashdot',
        label=r'$m = {}$'.format(m2))

# Perfect sinc repsonses
ax.plot(t, shifted_sinc, color='k', ls='dotted', alpha=0.75,
        label='$\mathrm{sinc} ( N \delta\omega t - \pi )$')

# Legend
ax.legend()

#--------------------------------#
#     Arrows Labelling Lines     #
#--------------------------------#
# # Arrow properties
# arrow_props = {'arrowstyle': '->', 'facecolor': 'k'}

# # m = 0 line
# ax.annotate(text=r'$m = 0$',
#             xytext=(0.4, 0.95), xy=(0.474, 0.769), 
#             arrowprops=arrow_props)

# # m = 1 line
# ax.annotate(text=r'$m = 1$',
#             xytext=(3.0, 0.7), xy=(1.719, 0.676),
#             arrowprops=arrow_props)

# # m = 2 line
# ax.annotate(text=r'$m = 2$',
#             xytext=(4.452, 0.450), xy=(4.024, 0.247),
#             arrowprops=arrow_props)

# # Sinc line
# ax.annotate(text=r'$\mathrm{sinc} \left( N \delta\omega t - m \right)$',
#             xytext=(2.7, 0.9), xy=(2.111, 0.821),
#             arrowprops=arrow_props)

#--------------------#
#     Axis Ticks     #
#--------------------#
ax.set_xticks(np.arange(0, max(t), 0.5*np.pi))
ax.set_xticks(np.arange(0.25*np.pi, max(t), 0.5*np.pi), minor=True)
ax.set_xticklabels(labels=[r'0', r'$\pi / 2$', r'$\pi$', r'$3 \pi / 2$',
                            r'$2 \pi$', r'$5 \pi / 2$', r'$3 \pi$'])

ax.set_yticks(np.arange(-0.2, 1.2, 0.2))
ax.set_yticks(np.arange(-0.1, 1.3, 0.2), minor=True)

#---------------------#
#     Axis limits     #
#---------------------#
ax.set_xlim(0, 3*np.pi)
ax.set_ylim(-0.25, 1.05)

#---------------------#
#     Axis Labels     #
#---------------------#
# Labels
ax.set_xlabel(r'$|\mathcal{E}_{d}| t$')
ax.set_ylabel(r'$\bar{A}(t)$')

#----------------------#
#     Figure stuff     #
#----------------------#
# Grid
ax.grid()

# fig.tight_layout(pad=0)
fig.savefig(filename_out)
fig.show()
