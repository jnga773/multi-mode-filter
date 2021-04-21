# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 19:10:04 2020

@author: Jacob
"""

import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.close('all')

import sys
# Append path to add Jacobs_Function
# sys.path.append("/home/jnga773/Documents/898-PhD/test/")
sys.path.append("E:/Google Drive/University/Physics/898 PhD/test")
from Jacobs_Functions import filename

#---------------------------------------#
#     Parameters from Fortran files     #
#---------------------------------------#
# Pull parameters from file
gamma, Omega, w0, kappa, epsilon, N, dw, phase, dt, t_max, tau1_max, tau2_max = \
    np.genfromtxt(fname=filename("parameters", "../data_files/"), delimiter="=", skip_header=1, usecols=(1,))
N = int(N)
phase = int(phase)
filters = (2*N) + 1

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

#---------------------#
#     QuTiP Stuff     #
#---------------------#
# Annihilation operators
# Mode frequency of modes a, b, and c
wa = w0 - dw
wb = w0
wc = w0 + dw
# Driving amplitude
Ea = np.sqrt(epsilon * gamma * kappa / filters) * np.exp(-1j * phase * np.pi / N)
Eb = np.sqrt(epsilon * gamma * kappa / filters)
Ec = np.sqrt(epsilon * gamma * kappa / filters) * np.exp(1j * phase * np.pi / N)
# Fock state truncation
Fock = 5

# Annihilation operators
psi0 = qt.tensor(qt.basis(Fock+1, 0), qt.basis(Fock+1, 0), qt.basis(Fock+1, 0), qt.basis(2, 1))

# Mode annihilation operators
a = qt.tensor(qt.destroy(Fock+1), qt.identity(Fock+1), qt.identity(Fock+1), qt.identity(2))
b = qt.tensor(qt.identity(Fock+1), qt.destroy(Fock+1), qt.identity(Fock+1), qt.identity(2))
c = qt.tensor(qt.identity(Fock+1), qt.identity(Fock+1), qt.destroy(Fock+1), qt.identity(2))
# Total cavity annihilation operator
A = a + b + c
# Atomic lowering operator
sm = qt.tensor(qt.identity(Fock+1), qt.identity(Fock+1), qt.identity(Fock+1), qt.sigmam())
sz = qt.tensor(qt.identity(Fock+1), qt.identity(Fock+1), qt.identity(Fock+1), qt.sigmaz())

# Hamiltonian
H = 0.5 * Omega * (sm + sm.dag())
H += (wa * a.dag() * a) + (wb * b.dag() * b) + (wc * c.dag() * c)
H += 0.5j * (Ea * a * sm.dag() + Eb * b * sm.dag() + Ec * c * sm.dag() - \
      np.conj(Ea) * a.dag() * sm - np.conj(Eb) * b.dag() * sm - np.conj(Ec) * c.dag() * sm)

# Collapse operators
c_ops = [np.sqrt(gamma * (1 - epsilon)) * sm,
         np.sqrt(kappa) * a,
         np.sqrt(epsilon * gamma / filters) * np.exp(-1j * phase * np.pi / N) * sm + np.sqrt(kappa) * a,
         np.sqrt(kappa) * b,
         np.sqrt(epsilon * gamma / filters) * sm + np.sqrt(kappa) * b,
         np.sqrt(kappa) * c,
         np.sqrt(epsilon * gamma / filters) * np.exp(1j * phase * np.pi / N) * sm + np.sqrt(kappa) * c]

# Evaluation operators
e_ops = [sz, A.dag() * A,
         a.dag() * b * c,
         a.dag() * b.dag() * c,
         a.dag() * b.dag() * c * b]

#--------------------------------#
#     Plot State Populations     #
#--------------------------------#
filename_states = filename("states", "../data_files/")
t = np.genfromtxt(filename_states, usecols=(0,))
popinv_f = np.genfromtxt(filename_states, usecols=(2,)) - np.genfromtxt(filename_states, usecols=(1,))
photon_f = np.genfromtxt(filename_states, usecols=(3,))
del filename_states

if N == 1:
    # Solve master equation
    data = qt.mesolve(H, psi0, t, c_ops, e_ops, progress_bar=True)
    popinv_q = data.expect[0]
    photon_q = data.expect[1]
    ataa = data.expect[2]
    atata = data.expect[3]

    # Plot
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=[12, 8])

    ax[0].plot(t, popinv_f, color='C0', ls='solid', label='Fortran')
    ax[0].plot(t, popinv_q, color='C1', ls='--', label='QuTiP')
    ax[0].set_ylabel(r'$ \langle \hat{\sigma}_{z} \rangle $', fontsize=15)

    ax[1].plot(t, photon_f, color='C0', ls='solid')
    ax[1].plot(t, photon_q, color='C1', ls='--')
    ax[1].set_ylabel(r'$ \langle \hat{a}^{\dagger} \hat{a} \rangle $', fontsize=15)
    ax[1].set_xlabel(r'$ \gamma \tau $', fontsize=15)
    ax[0].legend(loc='best', fontsize=15)

    ax[0].set_title(r'Fortran: $\left(\Omega = {}\gamma, \omega_0 = {}\gamma, \epsilon={}, \kappa = {}\gamma, N = {}, \delta\omega = {}\gamma \right)$'.format(Omega_str, w0_str, epsilon, kappa, N, dw), fontsize=15)
    fig.tight_layout()
    plt.show()

    print("Steady state photon number (Fortran) = {}".format(photon_f[-1]))
    print("Steady state photon number (QuTiP) = {}".format(photon_q[-1]))
    print(" ")
    print("<at_{{-1}} a_{{0}} a_{{1}}>_ss = {}".format(ataa[-1]))
    print("<at_{{-1}} at_{{0}} a_{{1}}>_ss = {}".format(atata[-1]))
