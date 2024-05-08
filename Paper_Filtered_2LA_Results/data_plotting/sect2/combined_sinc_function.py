#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 13:35:07 2023

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../paper_style.mplstyle')

plt.close('all')

# Filename for figure
filename_out = "../../images/sect2/fig1_combined_sinc_function.pdf"

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def sinc(t_in, phase):
    """
    Re-defining np.sinc function
    Parameters
    ----------
    t_in : array_like
        Times.
    Returns
    -------
    sinc : array_like
        sin(t - \phi) / (t - \phi)
    """
    from numpy import pi
    from numpy import sinc as npsinc
    sinc_f = npsinc((t_in - phase) / pi) / ((t_in - phase) / pi)
    return sinc_f

def theta(t_in):
    """
    Heaviside function
    """
    return np.heaviside(t_in, 0.5)

# Rectangle distribution
def rectangle(omega_in, a_in):
    # return (np.abs(omega_in)<=0.5*a_in).astype(float)
    rect = (0.5 * a_in) * (np.sign(a_in - omega_in) + np.sign(a_in + omega_in))
    return rect

def sinc_box(omega_in, a_in, phase_in):
    rect = (1 / a_in) * np.exp(1j * phase_in * omega_in / a_in)
    # Rectangle bit
    rect *= rectangle(omega_in, a_in)
    return rect

def sinc_fft(t_in, x_in):
    # Calculates the FFT of the input sinc function
    fft = np.fft.fft(x_in)
    norm = np.abs(fft[0])
    freq = np.fft.fftfreq(t_in.shape[0], t_in[1] - t_in[0])

    # Shift so arrays are numerically ordered
    fft = np.fft.fftshift(fft)
    freq = np.fft.fftshift(freq)
    # Power spectrum 
    fft = np.abs(fft) / norm
    # Multiply frequency by two
    freq *= 2.0

    # half way for checking
    just_before_zero = int(0.5 * len(x_in)) - 100
    if x_in[just_before_zero] != 0.0:
        fft = rectangle(freq, a)
        # fft = sinc_fft(freq, a, phase).real

    # Normalise to integral
    # normalise = np.sum(fft) * (freq[1] - freq[0])
    # fft *= (1 / normalise)

    # Normalise to max = 0.5
    fft = 0.5 * (fft / max(fft))
    
    return freq, fft

#-----------------------------------------------------------------------------#
#                                 PARAMETERS                                  #
#-----------------------------------------------------------------------------#
# Half-width of box
a = 1.0

# Time step spacing
dt = 0.01
# Time array
t = np.arange(-100, 100+dt, dt)

#-----------------------------------------------------------------------------#
#                             CALCULATE REPSONSES                             #
#-----------------------------------------------------------------------------#
# Time series: Sinc functions
time_series = [np.sinc(a * t),
               np.sinc(a * t) * theta(t),
               np.sinc(a * t - 1) * theta(t)]

# Frequency series: Fourier transforms
freq, fft0 = sinc_fft(t, time_series[0])
freq, fft1 = sinc_fft(t, time_series[1])
freq, fft2 = sinc_fft(t, time_series[2])
# Put into array
freq_series = [fft0,
               fft1,
               fft2]

#-----------------------------------------------------------------------------#
#                                 WRITE DATA                                  #
#-----------------------------------------------------------------------------#
# # Filenames to write data to
# filename_time = './data_time.txt'
# filename_freq = './data_freq.txt'

# # Write time data
# with open(filename_time, 'w') as file:
#     # Cycle through data indicies
#     for i in range(len(t)):
#         # String to write to file
#         str_write = ("{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t {:>20.15f} \n"
#                      ).format(t[i], time_series[0][i], time_series[1][i], time_series[2][i])
        
#         # Write to file
#         file.write(str_write)
        
# # Close file
# file.close()

# # Write frequency data
# with open(filename_freq, 'w') as file:
#     # Cycle through data indicies
#     for i in range(len(freq)):
#         # String to write to file
#         str_write = ("{:>20.15f} \t {:>20.15f} \t {:>20.15f} \t {:>20.15f} \n"
#                      ).format(t[i], freq_series[0][i], freq_series[1][i], freq_series[2][i])
        
#         # Write to file
#         file.write(str_write)
        
# # Close file
# file.close()

#-----------------------------------------------------------------------------#
#                                    PLOT                                     #
#-----------------------------------------------------------------------------#
# fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=False,
#                        figsize=[6, 7.5])
fig, ax = plt.subplots(nrows=3, ncols=2, sharex=False, sharey=False,
                       figsize=[6, 7])

for i in range(3):
    #--------------------------#
    #     Plot: Time Series    #
    #--------------------------#
    ax[i, 0].plot(t, time_series[i], color='k', ls='solid')

    #---------------------------------#
    #     Plot: Fourier Transform     #
    #---------------------------------#
    ax[i, 1].plot(freq, freq_series[i], color='k', ls='solid')

    #---------------------------------#
    #     Axis Ticks: Time Series     #
    #---------------------------------#
    # X-ticks: Major
    ax[i, 0].set_xticks(ticks=[-10, -5, 0, 5, 10])
    ax[i, 0].set_xticklabels(labels=[])
    # X-ticks: Minor
    ax[i, 0].set_xticks(ticks=[-7.5, -2.5, 0, 2.5, 7.5], minor=True)

    # Y-tick: Major
    ax[i, 0].set_yticks(ticks=[-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax[i, 0].set_yticks(ticks=[-0.1, 0.1, 0.3, 0.5, 0.7, 0.9],
                        minor=True)

    #---------------------------------------#
    #     Axis Ticks: Fourier Transform     #
    #---------------------------------------#
    # X-ticks: Major
    ax[i, 1].set_xticks(ticks=[-2, -1, 0, 1, 2])
    ax[i, 1].set_xticklabels(labels=[])
    # X-ticks: Minor
    ax[i, 1].set_xticks(ticks=[-1.5, -0.5, 0.5, 1.5],
                        minor=True)

    # Y-tick
    ax[i, 1].set_yticks(ticks=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
    ax[i, 1].set_yticks(ticks=[0.05, 0.15, 0.25, 0.35, 0.45],
                        minor=True)

    #---------------------#
    #     Axis Limits     #
    #---------------------#
    ax[i, 0].set_xlim(-10, 10)
    ax[i, 1].set_xlim(-2, 2)

    #----------------------#
    #     Figure Stuff     #
    #----------------------#
    # Grid
    ax[i, 0].grid()
    ax[i, 1].grid()

    #---------------------#
    #     Axis Labels     #
    #---------------------#
    ax[i, 1].set_ylabel('Fourier Transform')
    
#--------------------------#
#     Axis Tick Labels     #
#--------------------------#
ax[2, 0].set_xticklabels(labels=[r'-10$\pi$', r'-5$\pi$', r'0', r'5$\pi$', r'10$\pi$'],
                         minor=False)
ax[2, 1].set_xticklabels(labels=[r'$-2a$', r'$-a$', r'$0$', r'$a$', r'$2a$'],
                         minor=False)

#---------------------#
#     Axis Labels     #
#---------------------#
ax[0, 0].set_ylabel(r'$\mathrm{{sinc}} \left( a t \right)$')
ax[1, 0].set_ylabel(r'$\theta(t) \mathrm{{sinc}} \left( a t \right)$')
ax[2, 0].set_ylabel(r'$\theta(t) \mathrm{{sinc}} \left( a t - \pi \right)$')

ax[2, 0].set_xlabel(r'$a t$')
ax[2, 1].set_xlabel(r'$\omega$')

#---------------------#
#     Figure Text     #
#---------------------#
# Time series
x = -9.25
y = 0.9
ax[0, 0].text(x, y, s='(a)')
ax[1, 0].text(x, y, s='(c)')
ax[2, 0].text(x, y, s='(e)')

# Fourier transform
x = -1.85
y = 0.46
ax[0, 1].text(x, y, s='(b)')
ax[1, 1].text(x, y, s='(d)')
ax[2, 1].text(x, y, s='(f)')

#----------------------#
#     Figure Stuff     #
#----------------------#
# Figure stuff
# fig.tight_layout(pad=0)
fig.savefig(filename_out)
fig.show()
