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
from Jacobs_Functions import filename

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

#-----------------------------------------------------------------------------#
#                               FILENAME THINGS                               #
#-----------------------------------------------------------------------------#
# Read parameters
gamma, Omega, alpha, delta, xi, w0, kappa, dw, epsilon, N, phase = \
    np.genfromtxt(filename("parameters"), delimiter="=", skip_header=1, usecols=1)
xi = round(xi, 2)
N = int(N)
phase = int(phase)
Omega = round(Omega, 2)
w0 = round(w0, 2)


# Print out dressed state frequencies for these parameters
ew0, ewp, ewm = three_level_eig(Omega, alpha, delta, xi)
print(('w0 = {} \n'
       'wp = {} \n'
       'wm = {} \n').format(ew0, ewp, ewm))

# Read the data
tau = np.genfromtxt(filename("g2_corr"), usecols=0)
g2 = np.genfromtxt(filename("g2_corr"), usecols=1)

# Master equation to compare
# g2_me = np.genfromtxt(filename("me_g2_corr"), usecols=1)

# Unfiltered correlation
g2_atom = np.genfromtxt(filename("atom_g2_corr"), usecols=1)

#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
plt.figure(figsize=[10, 6])

plt.plot(tau, g2, color='C0', label='Filtered Correlation')
# plt.plot(tau, g2_me, color='C1', ls='--', label='Master Equation')

# Plot unfilterd correlation
plt.plot(tau, g2_atom, color='k', alpha=0.5, ls='dashed', label='Unfiltered Correlation')

# Plot dressed_state correlation
# plt.plot(tau, corr_dressed, ls='--', color='k', alpha=0.5, label='Dressed-State Approximation')

plt.xlabel(r'$\gamma \tau$', fontsize=15)
plt.ylabel(r'$g^{(2)}(\tau)$', fontsize=15)
plt.legend(loc='best', fontsize=15)
plt.title(r'$\left( \Omega = {}, \alpha = {}, \delta = {}, \xi = {} \right), \left( \kappa = {}, \omega_{{0}} = {}, N = {}, \delta\omega = {} \right)$'.format(Omega, alpha, delta, xi, kappa, w0 ,N, dw), fontsize=15)
plt.tight_layout()
plt.show()
