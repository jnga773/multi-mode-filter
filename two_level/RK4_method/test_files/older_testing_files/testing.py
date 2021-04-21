#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 11:26:27 2020

@author: jnga773
"""

import numpy as np
from python_analytic_moments import *

def asm_ss(gamma_in, Omega_in, j_in, kappa_in, w0_in, epsilon_in, N_in, dw_in):
    from python_analytic_moments import _wj, _Ej
    # parameters
    wj = _wj(w0_in, dw_in, j_in)
    Ej = _Ej(gamma_in, kappa_in, epsilon_in, N_in, j_in)
    alpha = -(0.5 * gamma_in + kappa_in + 1j * wj)
    # steady states
    sm = -1j * gamma_in * Omega_in / (2 * (Omega_in ** 2) + (gamma_in ** 2))
    sz = -(gamma_in ** 2) / (2 * (Omega_in ** 2) + (gamma_in ** 2))
    aj = -Ej * sm / (kappa_in + 1j * wj)
    # denominator
    denom = (2 * alpha) * (2 * (Omega_in ** 2) - (gamma_in * alpha) + 2 * (alpha ** 2))
    # calculate
    asm_out  = 2j * Omega_in * Ej * alpha * sm
    asm_out += (Omega_in ** 2) * Ej * (sz + 1)
    asm_out += -2j * gamma_in * Omega_in * alpha * aj
    asm_out *= (1 / denom)
    return asm_out

def asp_ss(gamma_in, Omega_in, j_in, kappa_in, w0_in, epsilon_in, N_in, dw_in):
    from python_analytic_moments import _wj, _Ej
    # parameters
    wj = _wj(w0_in, dw_in, j_in)
    Ej = _Ej(gamma_in, kappa_in, epsilon_in, N_in, j_in)
    alpha = -(0.5 * gamma_in + kappa_in + 1j * wj)
    # steady states
    sm = -1j * gamma_in * Omega_in / (2 * (Omega_in ** 2) + (gamma_in ** 2))
    sz = -(gamma_in ** 2) / (2 * (Omega_in ** 2) + (gamma_in ** 2))
    aj = -Ej * sm / (kappa_in + 1j * wj)
    # denominator
    denom = (2 * alpha) * (2 * (Omega_in ** 2) - (gamma_in * alpha) + 2 * (alpha ** 2))
    # calculate
    asp_out  = -2j * Omega_in * Ej * alpha * sm
    asp_out += Ej * ((Omega_in ** 2) + 2 * (alpha ** 2) - gamma_in * alpha) * (sz + 1)
    asp_out += 2j * gamma_in * Omega_in * alpha * aj
    asp_out *= (1 / denom)
    return asp_out

def asz_ss(gamma_in, Omega_in, j_in, kappa_in, w0_in, epsilon_in, N_in, dw_in):
    from python_analytic_moments import _wj, _Ej
    # parameters
    wj = _wj(w0_in, dw_in, j_in)
    Ej = _Ej(gamma_in, kappa_in, epsilon_in, N_in, j_in)
    alpha = -(0.5 * gamma_in + kappa_in + 1j * wj)
    # steady states
    sm = -1j * gamma_in * Omega_in / (2 * (Omega_in ** 2) + (gamma_in ** 2))
    sz = -(gamma_in ** 2) / (2 * (Omega_in ** 2) + (gamma_in ** 2))
    aj = -Ej * sm / (kappa_in + 1j * wj)
    # denominator
    denom = (2 * alpha) * (2 * (Omega_in ** 2) - (gamma_in * alpha) + 2 * (alpha ** 2))
    # calculate
    asz_out  = -4 * Ej * (alpha ** 2) * sm
    asz_out += 2j * Omega_in * Ej * alpha * (sz + 1)
    asz_out += 4 * gamma_in * (alpha ** 2) * aj
    asz_out *= (1 / denom)
    return asz_out

gamma = 1.0
Omega = 5.0 * np.pi
kappa = 2.0
dw = 1.0
N = 20
epsilon = 1.0
w0 = 0.0
j = 4

# matrix steady states
asm_mat, asp_mat, asz_mat = cavsig2_ss(gamma, Omega, j, kappa, w0, epsilon, N, dw, "a")
# Gross equation steady states
asm_eqn = asm_ss(gamma, Omega, j, kappa, w0, epsilon, N, dw)
asp_eqn = asp_ss(gamma, Omega, j, kappa, w0, epsilon, N, dw)
asz_eqn = asz_ss(gamma, Omega, j, kappa, w0, epsilon, N, dw)

print("Matrix inversion method:")
print(r"<a_{} \sigma_->_ss = ".format(j), asm_mat)
print(r"<a_{} \sigma_+>_ss = ".format(j), asp_mat)
print(r"<a_{} \sigma_z>_ss = ".format(j), asz_mat)
print("------------------------")
print("Gross Equation:")
print(r"<a_{} \sigma_->_ss = ".format(j), asm_eqn)
print(r"<a_{} \sigma_+>_ss = ".format(j), asp_eqn)
print(r"<a_{} \sigma_z>_ss = ".format(j), asz_eqn)