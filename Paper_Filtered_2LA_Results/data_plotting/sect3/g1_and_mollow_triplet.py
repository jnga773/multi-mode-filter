#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 12:44:15 2022

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../paper_style.mplstyle')

plt.close('all')

# Filename: g1 figure
# filename_g1 = "../../images/sect3/g1_low_high_drive.pdf"
# Filename: spectrum
filename_spec = "../../images/sect3/fig5_spectrum_low_high_drive.pdf"

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

def calc_g1(gamma_in, Omega_in, tau_in):
    """
    Calculate the first-order correlation function
        < \sigma_{+}(0) \sigma_{-}(\tau) >.
    """
    # Make the input parameters complex
    Omega_in += 0.0j
    gamma_in += 0.0j
    
    # Lil cheaty bits
    Y = np.sqrt(2) * Omega_in / gamma_in
    Y2 = Y ** 2
    delta = np.sqrt( ((0.25 * gamma_in) ** 2) - (Omega_in ** 2) )
    
    # First-order correlation function
    # Break it up into parts
    c1 = (0.25 * (Y2 / (1 + Y2))) * np.exp(-0.5 * gamma_in * tau_in)
    
    c2 = -(0.125 * (Y2 / ((1 + Y2) ** 2)))
    c2 *= (1 - Y2 + ((1 - 5 * Y2) * (0.25 * gamma_in) / delta))
    c2 *= np.exp(-((0.75 * gamma_in) - delta) * tau_in)
    
    c3 = -(0.125 * (Y2 / ((1 + Y2) ** 2)))
    c3 *= (1 - Y2 - ((1 - 5 * Y2) * (0.25 * gamma_in) / delta))
    c3 *= np.exp(-((0.75 * gamma_in) + delta) * tau_in)
    
    g1_out = c1 + c2 + c3
              
    # Normalise
    g1_out *= 1 / g1_out[0]
              
    return g1_out

#-----------------------------------------------------------------------------#
#                                  PARAMETERS                                 #
#-----------------------------------------------------------------------------#
# Atomic decay rate
gamma = 1.0
# Driving amplitude
Omega_lo = 0.1
Omega_hi = 7.5

# Time-step
dt = 0.001
# g1 - time
tau_max = 100.0
tau = np.arange(0, tau_max + dt, dt)

#-----------------------------------------------------------------------------#
#                   PLOT: First-Order Correlation Function                    #
#-----------------------------------------------------------------------------#
# Calculate the first-order correlation function
g1_lo = calc_g1(gamma, Omega_lo, tau)
g1_hi = calc_g1(gamma, Omega_hi, tau)

#--------------#
#     PLOT     #
#--------------#
# fig, ax = plt.subplots(num='First-Order Correlation Function',
#                        nrows=1, ncols=2, sharey=True)

# #-------------------------#
# #     Plot: Low drive     #
# #-------------------------#
# ax[0].plot(tau, g1_lo.real, color='C0', ls='solid')

# #--------------------------#
# #     Plot: High drive     #
# #--------------------------#
# ax[1].plot(tau, g1_hi.real, color='C0', ls='solid')

# #--------------------#
# #     Axis Ticks     #
# #--------------------#
# ax[0].set_xticks(np.arange(0.0, 12.0, 2.0))
# ax[0].set_xticks(np.arange(1.0, 12.0, 2.0), minor=True)

# ax[1].set_xticks(np.arange(0.0, 12.0, 2.0))
# ax[1].set_xticks(np.arange(1.0, 12.0, 2.0), minor=True)

# ax[0].set_yticks(np.arange(0.0, 1.2, 0.2))
# ax[0].set_yticks(np.arange(0.1, 1.2, 0.2), minor=True)

# #---------------------#
# #     Axis Limits     #
# #---------------------#
# ax[0].set_xlim(-0.5, 10.5)
# ax[1].set_xlim(-0.5, 10.5)

# ax[0].set_ylim(-0.05, 1.05)

# #---------------------#
# #     Axis Labels     #
# #---------------------#
# ax[0].set_xlabel(r'$\gamma \tau$')
# ax[1].set_xlabel(r'$\gamma \tau$')

# ax[0].set_ylabel(r'$g^{(1)}_{ss}(\tau)$')


# # fig.tight_layout()
# # fig.savefig(filename_g1)
# fig.show()

#-----------------------------------------------------------------------------#
#                                MOLLOW TRIPLET                               #
#-----------------------------------------------------------------------------#
# Calculate the power spectrum
spec_lo, wlist = spectrum(tau, g1_lo)
spec_hi, wlist = spectrum(tau, g1_hi)

#--------------#
#     PLOT     #
#--------------#
fig, ax = plt.subplots(num='Power Spectrum', nrows=1, ncols=2,
                       sharey=True)

#-------------------------#
#     Plot: Low drive     #
#-------------------------#
ax[0].plot(wlist, spec_lo, color='k', ls='solid')

#--------------------------#
#     Plot: High drive     #
#--------------------------#
ax[1].plot(wlist, spec_hi, color='k', ls='solid')

#--------------------#
#     Axis Ticks     #
#--------------------#
ax[0].set_xticks(np.arange(-10.0, 15.0, 5.0))
ax[0].set_xticks(np.arange(-7.5, 15.0, 5.0), minor=True)

ax[1].set_xticks(np.arange(-10.0, 15.0, 5.0))
ax[1].set_xticks(np.arange(-7.5, 15.0, 5.0), minor=True)

ax[0].set_yticks(np.arange(0.0, 1.2, 0.2))
ax[0].set_yticks(np.arange(0.1, 1.2, 0.2), minor=True)

#---------------------#
#     Axis Limits     #
#---------------------#
ax[0].set_xlim(-10, 10)
ax[1].set_xlim(-10, 10)

ax[0].set_ylim(-0.025, 1.025)

#---------------------#
#     Axis Labels     #
#---------------------#
ax[0].set_xlabel(r'$\left( \omega - \omega_{A} \right) / \gamma$')
ax[1].set_xlabel(r'$\left( \omega - \omega_{A} \right) / \gamma$')

ax[0].set_ylabel(r'$S_{\mathrm{inc}}(\omega)$ (a.u.)')

#---------------------#
#     Figure Text     #
#---------------------#
ax[0].text(x=-9.5, y=0.95, s='(a)')
ax[1].text(x=-9.5, y=0.95, s='(b)')

#----------------------#
#     Figure Stuff     #
#----------------------#
# fig.tight_layout(pad=0)
fig.savefig(filename_spec)
fig.show()