# -*- coding: utf-8 -*-
"""
Created on Fri May 22 16:05:45 2020

@author: Jacob
"""

import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
from scipy.signal import find_peaks
import sys
# Append path to add Jacobs_Function
# sys.path.append("/home/jnga773/Documents/898-PhD/test/")
sys.path.append("E:/Google Drive/University/Physics/898 PhD/test")
from Jacobs_Functions import filename, spectrum

plt.close('all')

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

"""
This program simulates a driven two-level atom  coupled into a single-mode
filter cavity using QuTiP.
"""

#---------------------------------------#
#     Parameters from Fortran files     #
#---------------------------------------#
# Pull parameters from file
gamma, Omega, w0, kappa, epsilon, dt, t_max, tau1_max, tau2_max = \
    np.genfromtxt(fname=filename("parameters", "../data_files/"), delimiter="=", skip_header=1, usecols=(1,))
Fock = 10

#--------------------#
#     Parameters     #
#--------------------#
# # Atom decay rate
# gamma = 1.0
# # Atom driving strength
# Omega = 5.0 * np.pi
# # Cavity decay rate
# kappa = 10.0
# # Cavity resonance frequency
# w0 = Omega
# # Fraction of fluorescence sent to cavity
# epsilon = 1.0
# # Fock basis truncation
# Fock = 10
# # Time step
# dt = 0.001
# # Max time to integrate state solver for
# t_max = 20.0
# tau1_max = 100.0
# tau2_max = 0.0

#---------------#
#     Stuff     #
#---------------#
# Initial number of photons in cavity
N_init = 0
# List of times
t = np.arange(0, t_max + dt, dt)
tau1 = np.arange(0, tau1_max + dt, dt)
tau2 = np.arange(0, tau2_max + dt, dt)

# Check if Omega/D0 is a multiple of pi
if np.mod(np.round(Omega / np.pi, 1), 1) == 0.0 and Omega != 0.0:
    Omega_str = "{} \pi".format(int(Omega / np.pi))
    if w0 == Omega or w0 == -Omega:
        w0_str = Omega_str
    elif w0 == 0.0:
        w0_str = "0.0"
    else:
        w0_str = str(round(w0, 2))
else:
    Omega_str = str(round(Omega, 2))
    w0_str = str(round(w0, 2))

#-------------------------#
#     QUTIP Operators     #
#-------------------------#
# Initial state vector
psi0 = qt.tensor(qt.basis(2, 1), qt.basis(Fock+1, N_init))

# Cavity annihilation operator
a = qt.tensor(qt.qeye(2), qt.destroy(Fock+1))
# Atom decay operator
sm = qt.tensor(qt.sigmam(), qt.qeye(Fock+1))
# sigmaz operator
sz = qt.tensor(qt.sigmaz(), qt.qeye(Fock+1))

# Hamiltonian
# Two-level atom
H = 0.5 * Omega * (sm + sm.dag())
# Cavity mode
H += w0 * a.dag() * a
# Cascade system coupling
H += 0.5j * np.sqrt(epsilon * gamma * kappa) * (sm.dag() * a - a.dag() * sm)

# Collapse operators
c_ops = [np.sqrt(gamma * (1.0 - epsilon)) * sm,
         np.sqrt(epsilon * gamma) * sm + np.sqrt(kappa) * a,
         np.sqrt(kappa) * a]

#----------------------#
#     State Solver     #
#----------------------#
# # Evaluation operators
# e_ops = [sz, a.dag() * a, a]

# # Solve master equation
# data = qt.mesolve(H, psi0, t, c_ops, e_ops, progress_bar=True)
# popinv = data.expect[0]
# photon = data.expect[1]

# print("QuTiP: Steady state photon number = {}".format(photon[-1]))

fig, ax = plt.subplots(2, 1, sharex=True, figsize=[12, 8])

ax[0].plot(t, popinv)
ax[0].set_ylabel(r'$ \langle \hat{\sigma}_{z} \rangle $', fontsize=15)

ax[1].plot(t, photon)
ax[1].set_ylabel(r'$ \langle \hat{A}^{\dagger} \hat{A} \rangle $', fontsize=15)
ax[1].set_xlabel(r'$ \gamma \tau $', fontsize=15)

ax[0].set_title(r'QuTiP $\left( \Omega = {}, \omega_0 = {}, \epsilon = {}, \kappa = {}, \mathrm{{truncation}} = {}~\mathrm{{ photons}} \right)$'.format(Omega_str, w0_str, epsilon, kappa, Fock), fontsize=15)
fig.tight_layout()
plt.show()

#---------------------------#
#     Operator Momemnts     #
#---------------------------#
# Evaluation operators
# e_ops = [a.dag() * a, a, a.dag(), a.dag() * sz]
e_ops = [sz, a.dag() * a,
         # First-order moments
         # sm, sm.dag(), sz, a, a.dag()]
         # Second-order moments
          # a * a, a.dag() * a.dag(),
          # a * sm, a * sm.dag(), a * sz,
          # a.dag() * sm, a.dag() * sm.dag(), a.dag() * sz]
         # Third-order moments
           # a.dag() * a * a, a.dag() * a.dag() * a,
           # a * a * sm, a * a * sm.dag(), a * a * sz,
           # a.dag() * a.dag() * sm, a.dag() * a.dag() * sm.dag(), a.dag() * a.dag() * sz,
           # a.dag() * a * sm, a.dag() * a * sm.dag(), a.dag() * a * sz]
         # Fourth-order moments
          a.dag() * a.dag() * a * a,
          a.dag() * a.dag() * a * sm, a.dag() * a.dag() * a * sm.dag(), a.dag() * a.dag() * a * sz,
          a.dag() * a * a * sm, a.dag() * a * a * sm.dag(), a.dag() * a * a * sz]


# Solve master equation
data = qt.mesolve(H, psi0, t, c_ops, e_ops, progress_bar=True)

popinv = data.expect[0]
photon = data.expect[1]

#---------------------#
# First-order moments #
#---------------------#
# _sm = data.expect[2]
# _sp = data.expect[3]
# _a = data.expect[4]
# _at = data.expect[5]
# print("--------------------------------")
# print("<sm>_ss = {}".format(_sm[-1]))
# print("<sp>_ss = {}".format(_sp[-1]))
# print("<sz>_ss = {}".format(popinv[-1]))
# print("--------------------------------")
# print("<a>_ss = {}".format(_a[-1]))
# print("<at>_ss = {}".format(_at[-1]))

#----------------------#
# Second-order moments #
#----------------------#
# aa = data.expect[2]
# atat = data.expect[3]
# asm = data.expect[4]
# asp = data.expect[5]
# asz = data.expect[6]
# atsm = data.expect[7]
# atsp = data.expect[8]
# atsz = data.expect[9]
# print("--------------------------------")
# print("<a a>_ss = {}".format(aa[-1]))
# print("<at at>_ss = {}".format(atat[-1]))
# print("--------------------------------")
# print("<a sm>_ss = {}".format(asm[-1]))
# print("<a sp>_ss = {}".format(asp[-1]))
# print("<a sz>_ss = {}".format(asz[-1]))
# print("--------------------------------")
# print("<at sm>_ss = {}".format(atsm[-1]))
# print("<at sp>_ss = {}".format(atsp[-1]))
# print("<at sz>_ss = {}".format(atsz[-1]))

#---------------------#
# Third-order moments #
#---------------------#
# ataa = data.expect[2]
# atata = data.expect[3]
# aasm = data.expect[4]
# aasp = data.expect[5]
# aasz = data.expect[6]
# atatsm = data.expect[7]
# atatsp = data.expect[8]
# atatsz = data.expect[9]
# atasm = data.expect[10]
# atasp = data.expect[11]
# atasz = data.expect[12]
# print("--------------------------------")
# print("<at a a>_ss = {}".format(ataa[-1]))
# print("<at at a>_ss = {}".format(atata[-1]))
# print("--------------------------------")
# print("<a a sm>_ss = {}".format(aasm[-1]))
# print("<a a sp>_ss = {}".format(aasp[-1]))
# print("<a a sz>_ss = {}".format(aasz[-1]))
# print("--------------------------------")
# print("<at at sm>_ss = {}".format(atatsm[-1]))
# print("<at at sp>_ss = {}".format(atatsp[-1]))
# print("<at at sz>_ss = {}".format(atatsz[-1]))
# print("--------------------------------")
# print("<at a sm>_ss = {}".format(atasm[-1]))
# print("<at a sp>_ss = {}".format(atasp[-1]))
# print("<at a sz>_ss = {}".format(atasz[-1]))

#----------------------#
# Fourth-order moments #
#----------------------#
atataa = data.expect[2]
atatasm = data.expect[3]
atatasp = data.expect[4]
atatasz = data.expect[5]
ataasm = data.expect[6]
ataasp = data.expect[7]
ataasz = data.expect[8]
print("--------------------------------")
print("<at at a a>_ss = {}".format(atataa[-1]))
print("--------------------------------")
print("<at at a sm>_ss = {}".format(atatasm[-1]))
print("<at at a sp>_ss = {}".format(atatasp[-1]))
print("<at at a sz>_ss = {}".format(atatasz[-1]))
print("--------------------------------")
print("<at a a sm>_ss = {}".format(ataasm[-1]))
print("<at a a sp>_ss = {}".format(ataasp[-1]))
print("<at a a sz>_ss = {}".format(ataasz[-1]))

#------------------------------------------#
#     First-Orger Correlation Function     #
#------------------------------------------#
if tau1_max > 0.0:
    # Solve first-order correlation function
    G1_filter = qt.correlation_2op_1t(H, None, tau1, c_ops, a_op=a.dag(), b_op=a)
    if photon[-1] != 0.0:
        g1_filter = G1_filter / photon[-1]
    else:
        g1_filter = G1_filter

    print("Mean Im[g1_cavity] = {}".format(np.mean(np.imag(g1_filter))))

    # Calculate spectrum from correlation data
    spec_filter, wlist_filter = spectrum(tau1, g1_filter)
    # Shift frequency by central resonance frequency of filter
    # wlist_multi = wlist_multi

    # Mollow triplet
    # Bare spectrum of atom (Mollow triplet)
    corr_atom = mollow_triplet(tau1, gamma, Omega)
    spec_atom, wlist_atom = spectrum(tau1, corr_atom)
    # Scale Mollow triplet so max is 1.0
    spec_atom = spec_atom / max(spec_atom)
    # Depending on where filter is centred, renormalise Mollow
    spec_atom_peaks = find_peaks(spec_atom, height=0.1)
    if np.abs(w0) <= 1.0:
        # Filter is centered on central peak
        spec_atom = max(spec_filter) * spec_atom / max(spec_atom)
    elif np.abs(w0) <= Omega+1 and np.abs(w0) >= Omega-1:
        # Filter is centered on side peak
        spec_atom = max(spec_filter) * spec_atom / spec_atom[spec_atom_peaks[0][0]]

    # Plot
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 8])
    # First-order correlation
    ax[0].plot(tau1, np.real(g1_filter), color='C0',
               label=r'$\mathrm{Re} \left( g^{(1)}(\tau) \right)$')
    ax[0].plot(tau1, np.imag(g1_filter), color='C1', lw=0.5, ls='--',
               label=r'$\mathrm{Im} \left( g^{(1)}(\tau) \right)$')
    ax[0].set_xlim(-1.0, 11.0)
    ax[0].set_xlabel(r'$\gamma \tau$', fontsize=15)
    ax[0].set_ylabel(r'$G^{(1)}(\tau)$', fontsize=15)
    ax[0].legend(loc='upper right', fontsize=15)

    # Plot spectrum
    ax[1].plot(wlist_filter, spec_filter, color='black', label='Cavity Spectrum')
    ax[1].plot(wlist_atom, spec_atom, color='k', ls='dotted', alpha=0.5, label='Scaled Mollow Triplet')
    ax[1].set_xlim(-int(1.5*Omega), int(1.5*Omega))
    ax[1].set_xlabel(r'$ \left( \omega - \omega_{d} \right) / \gamma $', fontsize=15)
    ax[1].set_ylabel(r'Power Spectrum', fontsize=15)
    ax[1].legend(loc='upper right', fontsize=15)

    fig.suptitle(r'QuTiP: $\left( \Omega={}, \omega_0={}, \kappa = {}, \epsilon={}, \mathrm{{truncation}} = {}~\mathrm{{photons}} \right)$'.format(Omega_str, w0_str, kappa, epsilon, Fock), fontsize=15)
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.show()

#-------------------------------------------#
#     Second-Order Correlation Function     #
#-------------------------------------------#
if tau2_max > 0.0:
    # Solve second-order correlation function
    G2 = qt.correlation_3op_1t(H, None, tau2, c_ops, a_op=a.dag(), b_op=a.dag() * a, c_op=a).real
    # G2 = qt.correlation_3op_1t(H, psi0, tau2, c_ops, a_op=a.dag(), b_op=a.dag() * a, c_op=a).real
    if photon[-1] != 0.0:
        g2 = G2 / (photon[-1] ** 2)
    else:
        g2 = G2

    print("Initial g2 value = {}".format(g2[0]))

    # e_decay = 0.5 * g2[0] * (1 + np.exp(-kappa * t))

    # Plot data
    plt.figure(figsize=[10, 6])
    plt.plot(tau2, g2)
    if w0 == 0.0:
        plt.plot(tau2, np.ones(len(tau2)), ls='--', color='k', alpha=0.5)
    elif np.abs(w0) == Omega:
        plt.plot(tau2, 1 - np.exp(-0.5*gamma*tau2), ls='--', color='k', alpha=0.5)
    # plt.plot(t, e_decay, color='C1', ls='--')
    plt.ylabel(r'$g^{(2)}(\tau)$', fontsize=15)
    plt.xlabel(r'$\kappa \tau$', fontsize=15)
    plt.title(r'QuTiP: $\left( \Omega={}, \omega_0={}, \kappa = {}, \epsilon={}, \mathrm{{truncation}} = {}~\mathrm{{photons}} \right)$'.format(Omega_str, w0_str, kappa, epsilon, Fock), fontsize=15)
    plt.tight_layout()
    plt.show()
