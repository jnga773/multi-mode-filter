# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 12:16:06 2020

@author: Jacob
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
import sys
# Append path to add Jacobs_Function
sys.path.append("/home/jnga773/Documents/898-PhD/test/")
# sys.path.append("D:/Google Drive/University/Physics/898 PhD/test")
from Jacobs_Functions import filename, spectrum

plt.close('all')

#-------------------#
#     Read Data     #
#-------------------#

# Filenames
# parent_dir = "./cluster/"
parent_dir = "../data_files/"
filename_param = filename("parameters", parent_dir)
filename_states = filename("states", parent_dir)
filename_mat_states = filename("matrix_states", parent_dir)
filename_g1 = filename("g1_corr", parent_dir)
filename_mat_g1 = filename("matrix_g1_corr", parent_dir)
filename_g2 = filename("g2_corr", parent_dir)
filename_mat_g2 = filename("matrix_g2_corr", parent_dir)

# Pull parameters from file
gamma, Omega, w0, kappa, epsilon, N, dw, phase, dt, t_max, tau1_max, tau2_max = \
    np.genfromtxt(fname=filename_param, delimiter="=", skip_header=1, usecols=(1,))
N = int(N)
phase = int(phase)
filters = (2*N) + 1
del(filename_param)

# Check if Omega/D0 is a multiple of pi
if np.mod(np.round(Omega / np.pi, 5), 1) == 0.0 and Omega != 0.0:
    Omega_str = "{} \pi".format(int(Omega / np.pi))
    if w0 == Omega:
        w0_str = Omega_str
    elif w0 == -Omega:
        w0_str = -Omega_str
    elif w0 == 0.0:
        w0_str = "0.0"
else:
    Omega_str = str(round(Omega,2))
    w0_str = str(round(w0, 2))

#--------------------------------#
#     Plot State Populations     #
#--------------------------------#
# Pull data from file
# time
t = np.genfromtxt(fname=filename_states, dtype='float', usecols=(0,))
# population inversions
popinv = np.genfromtxt(filename_states, usecols=(2,)) - np.genfromtxt(filename_states, usecols=(1,))
popinv_mat = np.genfromtxt(filename_mat_states, usecols=(2,)) - np.genfromtxt(filename_mat_states, usecols=(1,))
# mean photon number in cavities
photon = np.genfromtxt(fname=filename_states, dtype='float', usecols=(3,))
photon_mat = np.genfromtxt(filename_mat_states, usecols=(3,))
del filename_states, filename_mat_states

# Plot
fig, ax = plt.subplots(2, 1, sharex=True, figsize=[12, 8])

ax[0].plot(t, popinv, color='C0', ls='solid', label='Moment Equations')
ax[0].plot(t, popinv_mat, color='C1', ls='--', label='Matrix Equations')
ax[0].set_ylabel(r'$ \langle \hat{\sigma}_{z} \rangle $', fontsize=15)

ax[1].plot(t, photon, color='C0', ls='solid')
ax[1].plot(t, photon_mat, color='C1', ls='--')
ax[1].set_ylabel(r'$ \langle \hat{a}^{\dagger} \hat{a} \rangle $', fontsize=15)
ax[1].set_xlabel(r'$ \gamma \tau $', fontsize=15)
ax[0].legend(loc='best', fontsize=15)

ax[0].set_title(r'Fortran: $\left(\Omega = {}\gamma, \omega_0 = {}\gamma, \epsilon={}, \kappa = {}\gamma, N = {}, \delta\omega = {}\gamma \right)$'.format(Omega_str, w0_str, epsilon, kappa, N, dw), fontsize=15)
fig.tight_layout()
plt.show()

#-----------------------------------#
#     Plot G^{(1)} and Spectrum     #
#-----------------------------------#
if tau1_max > 0:
    # Pull data
    tau1 = np.genfromtxt(filename_g1, usecols=(0,))
    g1 = np.genfromtxt(filename_g1, usecols=(1,)) + 1j * np.genfromtxt(filename_g1, usecols=(2,))
    g1_mat = np.genfromtxt(filename_mat_g1, usecols=(1,)) + 1j * np.genfromtxt(filename_mat_g1, usecols=(2,))
    print("Mean Im[g1_moment] = {}".format(np.mean(g1.imag)))
    print("Mean Im[g1_matrix] = {}".format(np.mean(g1_mat.imag)))
    del filename_g1, filename_mat_g1

    def mollow_triplet(tau_in, gamma_in, Omega_in):
        """
        Expression from Howard vol. 1
        """
        # Turn inputs into complex values
        Omega_in = np.complex(Omega_in)
        gamma_in = np.complex(gamma_in)
        # input tau
        # Calculate first order correlation
        Yf = (np.sqrt(2.0) * Omega_in) / (gamma_in)
        df = np.sqrt(((0.25 * gamma_in) ** 2) - (Omega_in ** 2))
        g1_f = (0.25 * ((Yf ** 2) / (1 + (Yf ** 2)))) * np.exp(-0.5 * gamma_in * tau_in) - \
               (0.125 * ((Yf ** 2) / ((1 + (Yf ** 2)) ** 2))) * (1 - (Yf ** 2) + ((1 - 5 * (Yf ** 2)) * (0.25 * gamma_in) / df)) * np.exp(-((0.75 * gamma_in) - df) * tau_in) - \
               (0.125 * ((Yf ** 2) / ((1 + (Yf ** 2)) ** 2))) * (1 - (Yf ** 2) - ((1 - 5 * (Yf ** 2)) * (0.25 * gamma_in) / df)) * np.exp(-((0.75 * gamma_in) + df) * tau_in)
        print("Mean Im[g1_atom] = {}".format(np.mean(np.imag(g1_f))))
        return g1_f

    from scipy.signal import find_peaks

    spec, wlist = spectrum(tau1, g1, norm='peak')
    spec_mat, wlist = spectrum(tau1, g1_mat, norm='peak')
    spec_atom, wlist = spectrum(tau1, mollow_triplet(tau1, gamma, Omega), norm='peak')

    spec_atom_peaks = find_peaks(spec_atom, height=0.1)
    if np.abs(w0) <= 1.0:
        # Filter is centered on central peak
        spec_atom = max(spec) * spec_atom / max(spec_atom)
    elif np.abs(w0) <= Omega+1 and np.abs(w0) >= Omega-1:
        # Filter is centered on side peak
        spec_atom = max(spec) * spec_atom / spec_atom[spec_atom_peaks[0][0]]

    # Plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 8])

    # First-order correlation
    ax[0].plot(tau1, np.real(g1), color='C0', ls='solid', lw=2.0, alpha=1.0,
               label=r'Moment Equations (Real)')
    ax[0].plot(tau1, np.imag(g1)+1.0, color='C0', ls='dashed', lw=2.0, alpha=1.0,
               label=r'Moment Equations (Imag) + 1')
    ax[0].plot(tau1, np.real(g1_mat), color='C1', ls='solid', lw=0.75, alpha=1.0,
               label=r'Matrix Equations (Real)')
    ax[0].plot(tau1, np.imag(g1_mat)+1.0, color='C1', ls='dashed', lw=0.75, alpha=1.0,
               label=r'Matrix Equations (Imag) + 1')
    ax[0].set_xlim(-1.0, 11.0)
    ax[0].set_xlabel(r'$\gamma \tau$', fontsize=15)
    ax[0].set_ylabel(r'$g^{(1)}(\tau)$', fontsize=15)
    ax[0].legend(loc='upper right', fontsize=15)

    # Plot spectrum
    ax[1].plot(wlist, spec, color='C0', ls='solid', label='Moment Equations')
    ax[1].plot(wlist, spec_mat, color='C1', ls='dashed', label='Matrix Equations')
    ax[1].plot(wlist, spec_atom, color='k', ls='dotted', alpha=0.5, label='Scaled Mollow Triplet')
    ax[1].set_xlim(-int(1.5*Omega), int(1.5*Omega))
    ax[1].set_xlabel(r'$ \left( \omega - \omega_{d} \right) / \gamma $', fontsize=15)
    ax[1].set_ylabel(r'Power Spectrum', fontsize=15)
    ax[1].legend(loc='best', fontsize=15)

    fig.suptitle(r'Fortran: $\left(\Omega = {}\gamma, \omega_0 = {}\gamma, \epsilon={}, \kappa = {}\gamma, N = {}, \delta\omega = {}\gamma \right)$'.format(Omega_str, w0_str, epsilon, kappa, N, dw), fontsize=15)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.show()

#----------------------#
#     Plot G^{(2)}     #
#----------------------#
if tau2_max > 0:
    # Pull data
    tau2 = np.genfromtxt(filename_g2, usecols=(0,))
    g2 = np.genfromtxt(filename_g2, usecols=(1,))
    g2_mat = np.genfromtxt(filename_mat_g2, usecols=(1,))
    print("Initial g2 value (moment) = {}".format(g2[0]))
    print("Initial g2 value (matrix) = {}".format(g2_mat[0]))
    del filename_g2, filename_mat_g2

    # Plot
    plt.figure(figsize=[10, 6])
    plt.plot(tau2, g2, color='C0', ls='solid', label='Moment Equations')
    plt.plot(tau2, g2_mat, color='C1', ls='dashed', label='Matrix Solver')
    plt.xlabel(r'$\gamma \tau$', fontsize=15)
    plt.ylabel(r'$g^{(2)}(\tau)$', fontsize=15)
    plt.legend(loc='best', fontsize=15)
    plt.title(r'Fortran: $\left(\Omega = {}\gamma, \omega_0 = {}\gamma, \epsilon={}, \kappa = {}\gamma, N = {}, \delta\omega = {}\gamma \right)$'.format(Omega_str, w0_str, epsilon, kappa, N, dw), fontsize=15)
    plt.tight_layout()
    plt.show()
