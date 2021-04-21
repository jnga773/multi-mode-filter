# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

plt.close('all')

def filename(file_in):
    direct_f = "./data_files/"
    ext_f = ".txt"
    filename_out = direct_f + file_in + ext_f
    return filename_out

def mollow_triplet(Omega_in, gamma_in):
    """
    Expression from Howard vol. 1
    """
    # Turn inputs into complex values
    Omega_in = np.complex(Omega_in)
    gamma_in = np.complex(gamma_in)
    # input tau
    tau_f = np.arange(0, 500, 0.01)
    # Calculate first order correlation
    Yf = (np.sqrt(2.0) * Omega_in) / (gamma_in)
    df = np.sqrt(((0.25 * gamma_in) ** 2) - (Omega_in ** 2))
    g1_f = (0.25 * ((Yf ** 2) / (1 + (Yf ** 2)))) * np.exp(-0.5 * gamma_in * tau_f) - \
           (0.125 * ((Yf ** 2) / ((1 + (Yf ** 2)) ** 2))) * (1 - (Yf ** 2) + ((1 - 5 * (Yf ** 2)) * (0.25 * gamma_in) / df)) * np.exp(-((0.75 * gamma_in) - df) * tau_f) - \
           (0.125 * ((Yf ** 2) / ((1 + (Yf ** 2)) ** 2))) * (1 - (Yf ** 2) - ((1 - 5 * (Yf ** 2)) * (0.25 * gamma_in) / df)) * np.exp(-((0.75 * gamma_in) + df) * tau_f)
    # Fourier transform for power spectrum
    fft = np.fft.fft(g1_f)#, norm='ortho')
    fft = np.fft.fftshift(fft)
    freq = np.fft.fftfreq(tau_f.shape[0], tau_f[1]-tau_f[0])
    freq = np.fft.fftshift(freq)
    # Remove zero frequency term
    indices = np.where(freq != 0.0)[0]
    spec_output = fft[indices]
    # Divide by 2 \pi
#    spec_output = (np.abs(spec_output) ** 2) / (2 * np.pi)
    spec_output = spec_output.real
    # take away non zero tails
    spec_output = spec_output - np.mean(spec_output[0])
    # The spectrum is given by the absolute value of the FFT squared.
#    spec = np.abs(spec) ** 2
    wlist_output = freq[indices] # wlist is in terms of FFT frequencies
    wlist_output = np.array(wlist_output) * 2 * np.pi
    return spec_output, wlist_output

#------------------------------------------------------------------------------#
#                                FILENAME THINGS                               #
#------------------------------------------------------------------------------#
# Read parameters
gamma, Omega, kappa, dw, epsilon, N, phase = \
    np.genfromtxt(filename("parameters"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
Omega = round(Omega, 2)

# Pull data from file
# Central frequencies
w0_list = np.genfromtxt(filename("photon"), usecols=0)
# Steady state photon number
photon = np.genfromtxt(filename("photon"), usecols=1)

# Bare spectrum of atom (Mollow triplet)
spec_atom, wlist_atom = mollow_triplet(Omega, gamma)
# Mask such that only values around D0al are plotted
spec = spec_atom[np.abs(wlist_atom) < max(w0_list)]
wlist = wlist_atom[np.abs(wlist_atom) < max(w0_list)]
del spec_atom, wlist_atom

spec = max(photon) * spec / max(spec)

#------------------------------------------------------------------------------#
#                            PLOT PHOTON NUMBER SCAN                           #
#------------------------------------------------------------------------------#
plt.figure(figsize=[12, 8])

# Plot photon number scan
plt.plot(w0_list, photon)
# Plot bare spectrum of atom
plt.plot(wlist, spec, color='k', ls='dotted', alpha=0.5)

plt.ylabel(r'$ \langle A^{\dagger} A \rangle_{ss} $', fontsize=15)
plt.xlabel(r'$ \omega_{0} / \gamma $', fontsize=15)
plt.title(r'$ \Omega = %s, \mathrm{phase} = %s\pi, \kappa = %s, \delta\omega = %s, N = %s $'%(Omega, phase, kappa, dw, N), fontsize=15)

plt.tight_layout()
plt.show()
