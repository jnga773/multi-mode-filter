# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
import sys
# Append path to add Jacobs_Function
sys.path.append("/home/jnga773/Documents/898-PhD/test/")
# sys.path.append("D:/Jacob/Google Drive/University/Physics/898 PhD/test")
from Jacobs_Functions import filename, spectrum

plt.close('all')

# Three-level Hamiltonian eigenvalues
def three_level_eig(Omega_in, alpha_in, delta_in, xi_in, vals_or_vecs='vals'):
    """
    Calculates the dressed/eigenvalues of the Hamiltonian in matrix form.
    """
    from numpy.linalg import eig
    from numpy import matrix
    # Set up matrix
    H = matrix([[0, 0.5 * Omega_in, 0],
                [0.5 * Omega_in, -(0.5 * alpha_in) - delta_in, 0.5 * xi_in * Omega_in],
                [0, 0.5 * xi_in * Omega_in, -2 * delta_in]])
    # Calculate eigenvalues
    # eigvals_out = eigvals(H)
    eigvals_out, eigvecs_out = eig(H)

    # Get the indicies that would sort them from big to small
    sorted_indices = eigvals_out.argsort()[::-1]
    # Return the sorted eigenvalues and eigenvectors
    eigvals_out = eigvals_out[sorted_indices]
    eigvecs_out = eigvecs_out[:, sorted_indices]

    # This sorting will be in the order [\omega_{+}, \omega_{0}, \omega_{-}].
    # To match it with the analytic calculations I've done, let's change it to
    # the the order of [\omega_{0}, \omega_{+}, \omega_{-}]
    sorted_indices = [1, 0, 2]
    eigvals_out = eigvals_out[sorted_indices]
    eigvecs_out = eigvecs_out[:, sorted_indices]

    # Return the output depending on vals_or_vecs
    if vals_or_vecs == 'vals':
        return eigvals_out
    elif vals_or_vecs == 'vecs':
        return eigvecs_out
    elif vals_or_vecs == 'both':
        return eigvals_out, eigvecs_out

def Gamma_values(Omega_in, alpha_in, delta_in, xi_in, Gamma_in, diag_or_off="off"):
    """
    Calculates the dressed state Gamma values for each transition
    """
    from numpy.linalg import inv
    from numpy import matrix, array
    eigvals, eigvecs = three_level_eig(Omega_in, alpha_in, delta_in, xi_in, "both")
    
    # Invert eigvec matrix
    Sinv = inv(eigvecs)
    
    # Grab each element from this matrix
    g0 = Sinv[0, 0]; gp = Sinv[1, 0]; gm = Sinv[2, 0]
    e0 = Sinv[0, 1]; ep = Sinv[1, 1]; em = Sinv[2, 1]
    f0 = Sinv[0, 2]; fp = Sinv[1, 2]; fm = Sinv[2, 2]
    
    # Nine matrix elements
    a1 = (g0 * e0) + xi_in * (e0 * f0)
    a2 = (g0 * ep) + xi_in * (e0 * fp)
    a3 = (g0 * em) + xi_in * (e0 * fm)
    a4 = (gp * e0) + xi_in * (ep * f0)
    a5 = (gp * ep) + xi_in * (ep * fp)
    a6 = (gp * em) + xi_in * (ep * fm)
    a7 = (gm * e0) + xi_in * (em * f0)
    a8 = (gm * ep) + xi_in * (em * fp)
    a9 = (gm * em) + xi_in * (em * fm)
    
    
    if diag_or_off == "off":
        # Return the Gamma values for the off-diagonal operators/density matrix components
        # Different diagonal Gamma rates
        Gpz = Gamma_in * ((a1 ** 2) + (a2 ** 2) + (a4 ** 2) + (a5 ** 2) +
                          (a7 ** 2) + (a8 ** 2) - (a1 * a5))
        Gmz = Gamma_in * ((a1 ** 2) + (a3 ** 2) + (a4 ** 2) + (a6 ** 2) +
                          (a7 ** 2) + (a9 ** 2) - (a1 * a9))
        Gpm = Gamma_in * ((a2 ** 2) + (a3 ** 2) + (a5 ** 2) + (a6 ** 2) +
                          (a8 ** 2) + (a9 ** 2) - (a5 * a9))
    
        return Gpz, Gmz, Gpm

    elif diag_or_off == "diag":
        # Return the matrix and B vector for the coupled systems
        M = matrix([[-Gamma_in * ((a4 ** 2) + (a6 ** 2) + (a7 ** 2)), Gamma_in * ((a2 ** 2) - (a6 ** 2))],
                   [Gamma_in * ((a4 ** 2) - (a7 ** 2)), -Gamma_in * ((a2 ** 2) + (a7 ** 2) + (a8 ** 2))]])
        B = array([Gamma_in * (a6 ** 2), Gamma_in * (a7 ** 2)])
        
        return M, B
    
def g1_dressed(Omega_in, alpha_in, delta_in, xi_in, Gamma_in, tau_in,
               ket_in, bra_in, diag_or_off="off"):
    """
    Calculates the first-order correlation function for two transitions.
    """
    from numpy import exp, matrix
    # Get eigenvalues
    eigenvalues = three_level_eig(Omega_in, alpha_in, delta_in, xi_in, 'vals')
    
    # Get Gamma values 
    Gpz, Gmz, Gpm = Gamma_values(Omega_in, alpha_in, delta_in, xi_in, Gamma_in, "off")
    # Turn it into a matrix
    Gamma_mat = matrix([[0.0, Gpz, Gmz],
                        [Gpz, 0.0, Gpm],
                        [Gmz, Gpm, 0.0]])

    # List of possible ket/bra inputs
    dressed_state_list = ["0", "+", "-"]
    
    if ket_in != bra_in:
        for ket_index, ket in enumerate(dressed_state_list):
            for bra_index, bra in enumerate(dressed_state_list):
                if ket == ket_in and bra == bra_in:
                    # print("| {} >  -->  | {} >".format(ket, bra))
                    # print("Gamma_pz = {}".format(Gpz))
                    # print("Gamma    = {}".format(Gamma_mat[ket_index, bra_index]))
                    
                    # Grab Gamma value
                    Gamma_ij = Gamma_mat[ket_index, bra_index] 
                    # Grab eigenvalues
                    wi = eigenvalues[ket_index]
                    wj = eigenvalues[bra_index]
                    
                    # Calculate first-order correlation function
                    g1 = exp((-(0.5 * Gamma_ij) + 1j * (wi - wj)) * tau_in)
    else:
        from sys import exit as sysexit
        print("I haven't figured out the |i> --> |i> g1 yet")
        sysexit(69)
    
    return g1

    return spec_out
#------------------------------------------------------------------------------#
#                                FILENAME THINGS                               #
#------------------------------------------------------------------------------#
# Read parameters
Gamma, Omega, alpha, delta, xi, w0, kappa, dw, epsilon, N, phase = \
    np.genfromtxt(filename("parameters"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
Omega = round(Omega, 2)
w0 = round(w0, 2)
xi = round(xi, 2)

# Print out dressed state frequencies for these parameters
ew0, ewp, ewm = three_level_eig(Omega, alpha, delta, xi)
print(('w0 = {} \n'
       'wp = {} \n'
       'wm = {} \n').format(ew0, ewp, ewm))

# ket to bra transition
ket = "+"
bra = "-"

# Read data
tau = np.genfromtxt(filename("g1_corr"), usecols=0)
g1 = np.genfromtxt(filename("g1_corr"), usecols=1) + 1j * \
     np.genfromtxt(filename("g1_corr"), usecols=2)

# Atomic/unfiltered correlation
g1_approx = g1_dressed(Omega, alpha, delta, xi, Gamma, tau, ket, bra, diag_or_off="off")

# Calculate Fourier transform
spec, wlist = spectrum(tau, g1, norm='peak')
spec_dressed, wlist = spectrum(tau, g1_approx, norm='peak')

# Renormalise dressed state spectrum
spec_dressed *= max(spec) / max(spec_dressed)

#-----------------------------------------------------------------------------#
#                               PLOT G^{(1)}                                  #
#-----------------------------------------------------------------------------#
# fig, ax = plt.subplots(2, 1, sharex=True, figsize=[8, 6])

# # Real part
# ax[0].plot(tau, g1.real, color='C0', ls='solid', lw=2.0, label='Filtered Spectrum')
# ax[0].plot(tau, g1_me.real, color='C1', ls='dashed', lw=1.0, label='Full Spectrum')
# ax[0].set_ylabel(r'Real Part', fontsize=15)
# ax[0].set_xlim(-0.2, 5.2)
# ax[0].legend(loc='best', fontsize=15)
# ax[0].set_title(r'$G^{{(1)}}(\tau)$ with $\left( \Omega = {} \gamma, \alpha = {} \gamma, \delta = {} \gamma, \xi = {} \right)$'.format(Omega, alpha, delta, xi), fontsize=15)

# ax[1].plot(tau, g1.imag, color='C0', ls='solid', lw=2.0)
# ax[1].plot(tau, g1_me.imag, color='k', ls='dashed', lw=1.0, alpha=0.5)
# ax[1].set_ylabel(r'Imaginary Part', fontsize=15)
# ax[1].set_xlabel(r'$\gamma t$', fontsize=15)
# ax[1].set_title(r'$N = {}, \kappa = {} \gamma, \omega_{{0}} = {} \gamma, \delta\omega = {} \gamma, \epsilon = {}$'.format(N, kappa, w0, dw, epsilon), fontsize=15)


# fig.tight_layout()
# fig.show()

#-----------------------------------------------------------------------------#
#                               PLOT SPECTRUM                                 #
#-----------------------------------------------------------------------------#
plt.figure(figsize=[8, 8])

plt.plot(wlist, spec, color='C0', ls='solid', lw=2.0, label='Filtered Spectrum')
plt.plot(wlist, spec_dressed, color='k', ls='dashed', lw=1.0, alpha=0.5, label='Dressed State Correlation')

# Set lims for three-level or effective two-level
if 2 * delta == -alpha:
    plt.xlim(-1.5 * Omega, 1.5 * Omega)
else:
    plt.xlim(-abs(alpha), abs(alpha))
plt.xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=15)
plt.ylabel(r'Power Spectrum (a.u.)', fontsize=15)

plt.legend(loc='best', fontsize=15)
plt.title(r'$\Omega = {} \gamma, \alpha = {} \gamma, \delta = {} \gamma, \xi = {}, N = {}, \delta\omega = {} \gamma, \kappa = {} \gamma, \omega_{{0}} = {} \gamma$'.format(Omega, alpha, delta, xi, N, dw, kappa, w0))

plt.tight_layout()
plt.show()

#-----------------------------------------------------------------------------#
#                          PLOT G^{(1)} AND SPECTRUM                          #
#-----------------------------------------------------------------------------#
# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 8])

# # First-order correlation
# ax[0].plot(tau, g1.real, color='C0', ls='solid', lw=2.0, label='Real')
# ax[0].plot(tau, g1.imag, color='C1', ls='dashed', lw=1.0, label='Imaginary')
# ax[0].set_xlim(-0.1, 10.1)

# ax[0].set_xlabel(r'$\gamma \tau$', fontsize=12)
# ax[0].set_ylabel(r'$G^{(1)}(\tau)$', fontsize=12)
# ax[0].set_title(r'$\Omega = {} \gamma, \alpha = {} \gamma, \delta = {} \gamma, \xi = {}$'.format(Omega, alpha, delta, xi), fontsize=12)
# ax[0].legend(loc='best', fontsize=12)

# # Spectra
# ax[1].plot(wlist, spec, color='C0', ls='solid', lw=2.0, label='Filtered')
# ax[1].plot(wlist, spec_me, color='k', ls='dashed', lw=1.0, alpha=0.5, label='Unfiltered')

# # Set lims for three-level or effective two-level
# if 2 * delta == -alpha:
#     ax[1].set_xlim(-1.5 * Omega, 1.5 * Omega)
# else:
#     ax[1].set_xlim(-abs(alpha), abs(alpha))
# ax[1].set_xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=12)
# ax[1].set_ylabel('Power Spectrum (a.u.)', fontsize=12)
# ax[1].set_title(r'$N = {}, \delta\omega = {} \gamma, \kappa = {} \gamma, \omega_{{0}} = {} \gamma$'.format(N, dw, kappa, w0), fontsize=12)
# ax[1].legend(loc='best', fontsize=12)

# fig.tight_layout()
# fig.show()

#-----------------------------------------------------------------------------#
#                        SPECTRUM AND ZOOMED SPECTRUM                         #
#-----------------------------------------------------------------------------#
# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 8])

# # Full-picture spectra
# ax[0].plot(wlist, spec, color='C0', ls='solid', lw=2.0, label='Filtered')
# ax[0].plot(wlist, spec_me, color='k', ls='dashed', lw=1.0, alpha=0.5, label='Unfiltered')

# ax[0].set_xlim(-abs(alpha), abs(alpha))
# ax[0].set_xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=12)
# ax[0].set_ylabel('Power Spectrum (a.u.)', fontsize=12)
# ax[0].set_title(r'$N = {}, \delta\omega = {} \gamma, \kappa = {} \gamma, \omega_{{0}} = {} \gamma$'.format(N, dw, kappa, w0), fontsize=12)
# ax[0].legend(loc='best', fontsize=12)


# # Zoomed-in spectra
# ax[1].plot(wlist, spec, color='C0', ls='solid', lw=2.0)
# ax[1].plot(wlist, spec_me, color='k', ls='dashed', lw=1.0, alpha=0.5)

# # Set lims for three-level or effective two-level
# if 2 * delta == -alpha:
#     ax[1].set_xlim(-1.5 * Omega, 1.5 * Omega)
# else:
#     ax[1].set_xlim(-abs(alpha), abs(alpha))
# ax[1].set_ylim(-0.001, 0.1)
# ax[1].set_xlabel(r'$\left( \omega - \omega_{d} \right) / \gamma$', fontsize=12)
# ax[1].set_ylabel('Power Spectrum (a.u.)', fontsize=12)
# ax[1].set_title(r'The same but zoomed in a bit', fontsize=12)


# fig.tight_layout()
# fig.show()
