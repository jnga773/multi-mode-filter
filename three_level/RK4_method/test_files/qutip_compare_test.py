#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 12:50:19 2017

@author: jacob
"""

import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.close('all')

def filename(filename_in, parent_in="./data_files/", ext_in="txt"):
    filename_out = parent_in + filename_in + "." + ext_in
    return filename_out

#-----------------------------------------------------------------------------#
#                                FILENAME THINGS                              #
#-----------------------------------------------------------------------------#
# Read parameters
gamma, Omega, alpha, delta, xi, w0, kappa, dw, epsilon, N, phase, dt, t_max, tau1_max, tau2_max = \
    np.genfromtxt(filename("parameters", parent_in="../data_files/"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
filters = (2*N) + 1

# Read data
t = np.genfromtxt(filename("states", parent_in="../data_files/"), usecols=0)
gg = np.genfromtxt(filename("states", parent_in="../data_files/"), usecols=1)
ee = np.genfromtxt(filename("states", parent_in="../data_files/"), usecols=2)
ff = np.genfromtxt(filename("states", parent_in="../data_files/"), usecols=3)
photon = np.genfromtxt(filename("states", parent_in="../data_files/"), usecols=4)
# # Tried and tested master equation
# gg_me = np.genfromtxt(filename("me_states"), usecols=1)
# ee_me = np.genfromtxt(filename("me_states"), usecols=2)
# ff_me = np.genfromtxt(filename("me_states"), usecols=3)
# photon_me = np.genfromtxt(filename("me_states"), usecols=4)

# print("gg_ss = {:.5f}, gg_me_ss = {:.5f}".format(gg[-1], gg_me[-1]))
# print("ee_ss = {:.5f}, ff_me_ss = {:.5f}".format(ee[-1], ee_me[-1]))
# print("ff_ss = {:.5f}, ee_me_ss = {:.5f}".format(ff[-1], ff_me[-1]))
# print("photon_ss = {:.5f}, photon_me_ss = {:.5f}".format(photon[-1], photon_me[-1]))

# Operator moments
moment = np.genfromtxt(filename("operators", parent_in="../data_files/"), usecols=1) + 1j * \
          np.genfromtxt(filename("operators", parent_in="../data_files/"), usecols=2)
# moment_me = np.genfromtxt(filename("me_operators"), usecols=1) + 1j * \
#             np.genfromtxt(filename("me_operators"), usecols=2)

#-----------------------------------------------------------------------------#
#                                 QuTiP STUFF                                 #
#-----------------------------------------------------------------------------#
Fock = 10

# Three-level state basis |g>, |e>, |f>
g, e, f = [qt.basis(3, 0), qt.basis(3, 1), qt.basis(3, 2)]

# Atomic operators
sgg = g * g.dag()
sgg = qt.tensor(sgg, qt.qeye(Fock+1))
sge = g * e.dag()
sge = qt.tensor(sge, qt.qeye(Fock+1))
seg = e * g.dag()
seg = qt.tensor(seg, qt.qeye(Fock+1))
see = e * e.dag()
see = qt.tensor(see, qt.qeye(Fock+1))
sef = e * f.dag()
sef = qt.tensor(sef, qt.qeye(Fock+1))
sfe = f * e.dag()
sfe = qt.tensor(sfe, qt.qeye(Fock+1))
sgf = g * f.dag()
sgf = qt.tensor(sgf, qt.qeye(Fock+1))
sfg = f * g.dag()
sfg = qt.tensor(sfg, qt.qeye(Fock+1))
sff = f * f.dag()
sff = qt.tensor(sff, qt.qeye(Fock+1))
# Lowering and raising operators
sm = sge + xi * sef

# Cavity annihilation operator
a = qt.tensor(qt.qeye(3), qt.destroy(Fock+1))

# Hamilton
# Atom
H = -(0.5*alpha + delta) * see - 2.0*delta * sff + 0.5*Omega * (sm + sm.dag())
# Cavity
H += w0 * a.dag() * a
# Cascade system coupling
H += 0.5j * np.sqrt(epsilon * gamma * kappa) * (sm.dag() * a - a.dag() * sm)

# # Initial state
# # psi0 = qt.tensor(g, qt.basis(Fock+1, N_init))
# psi = []
# for i in range(Fock+1):
#     psi.append([])
#     if i != 0:
#         temp = i / (np.sqrt(i))
#     else:
#         temp = 0
#     psi[i].append(temp)
# psi = qt.Qobj(psi)
# psi0 = psi / psi.norm()

# psi0 = qt.tensor(e, psi0)
# psi0 = psi0 / psi0.norm()

psi0 = qt.tensor(g, qt.basis(Fock+1, 0))

# Collapse operators
c_ops = [np.sqrt(gamma * (1.0 - epsilon)) * sm,
         np.sqrt(epsilon * gamma) * sm + np.sqrt(kappa) * a,
         np.sqrt(kappa) * a]


# the moment operator
moment_operator = a.dag() * a.dag() * a * a
# Evaluation operators
e_ops = [sgg, see, sff, a.dag() * a, moment_operator]


# Solve master equation
data = qt.mesolve(H, psi0, t, c_ops, e_ops, progress_bar=True)

gg_me = data.expect[0]
ee_me = data.expect[1]
ff_me = data.expect[2]
photon_me = data.expect[3]
moment_me = data.expect[4]
# moment_me = moment_me / moment_me[0]

Omega = round(Omega, 3)
xi = round(xi, 3)
w0 = round(w0, 3)

#-----------------------------------------------------------------------------#
#                               PLOT STATES                                   #
#-----------------------------------------------------------------------------#

# fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=[12, 10])

# # < \sigma^{gg}(t) >
# ax[0].plot(t, gg, color='C0', ls='solid', lw=2.0, label='New Equations')
# ax[0].plot(t, gg_me, color='C1', ls='dashed', lw=1.0, label='9x9 Master Equation')
# ax[0].set_ylim(-0.05, 1.05)
# ax[0].set_ylabel(r'$\langle \hat{\sigma}^{gg}(t) \rangle$', fontsize=15)
# ax[0].set_title(r'Three-Level Atom $\left( \Omega = {} \gamma, \alpha = {} \gamma, \delta = {} \gamma, \xi = {} \right)$'.format(round(Omega, 2), alpha, delta, xi), fontsize=15)

# # < \sigma^{ee}(t) >
# ax[1].plot(t, ee, color='C0', ls='solid', lw=2.0)
# ax[1].plot(t, ee_me, color='C1', ls='dashed', lw=1.0)
# ax[1].set_ylim(-0.05, 1.05)
# ax[1].set_ylabel(r'$\langle \hat{\sigma}^{ee}(t) \rangle$', fontsize=15)

# # < \sigma^{ff}(t) >
# ax[2].plot(t, ff, color='C0', ls='solid', lw=2.0)
# ax[2].plot(t, ff_me, color='C1', ls='dashed', lw=1.0)
# ax[2].set_ylim(-0.05, 1.05)
# ax[2].set_ylabel(r'$\langle \hat{\sigma}^{ff}(t) \rangle$', fontsize=15)

# # < A^{\dagger} A(t) >
# ax[3].plot(t, photon, color='C0', ls='solid', lw=2.0)
# ax[3].plot(t, photon_me, color='C1', ls='dashed', lw=1.0)
# ax[3].set_ylabel(r'$\langle \hat{A}^{\dagger} \hat{A}(t) \rangle$', fontsize=15)
# ax[3].set_title(r'$N = {}, \kappa = {} \gamma, \omega_{{0}} = {} \gamma, \delta\omega = {} \gamma, \epsilon = {}$'.format(N, kappa, w0, dw, epsilon), fontsize=15)

# # Labels
# ax[0].legend(loc='best', fontsize=15)
# ax[3].set_xlabel(r'$\gamma t$', fontsize=15)

# fig.tight_layout()
# fig.show()

#-----------------------------------------------------------------------------#
#                              PLOT MOMENTS                                   #
#-----------------------------------------------------------------------------#

fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=[8, 6])

# Real part
ax[0].plot(t, moment.real, color='C0', ls='solid', lw=2.0, label='Fortran')
ax[0].plot(t, moment_me.real, color='C1', ls='dashed', lw=1.0, label='QuTiP')
ax[0].set_ylabel(r'Real Part', fontsize=15)
ax[0].legend(loc='best', fontsize=15)
ax[0].set_title(r'$\Gamma = {}, \Omega = {}, \alpha = {}, \delta = {}, \xi = {}$'.format(gamma, Omega, alpha, delta, xi), fontsize=15)


ax[1].plot(t, moment.imag, color='C0', ls='solid', lw=2.0)
ax[1].plot(t, moment_me.imag, color='C1', ls='dashed', lw=1.0)
ax[1].set_title(r'$N = {}, \kappa = {}, \omega_{{0}} = {}, \delta\omega = {}, \epsilon = {}$'.format(N, kappa, w0, dw, epsilon), fontsize=15)
ax[1].set_ylabel(r'Imaginary Part', fontsize=15)
ax[1].set_xlabel(r'$\gamma t$', fontsize=15)

fig.tight_layout()
fig.show()
