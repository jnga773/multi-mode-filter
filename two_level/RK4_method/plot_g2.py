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

def g2_analytic(tau_in, Omega_in, gamma_in):
    """
    Analytic expression for the normalised second-order correlation function
    from Howard 1
    """
    from numpy import exp, sqrt, cosh, sinh, real
    Omega_in = complex(Omega_in)
    gamma_in = complex(gamma_in)
    d_f = sqrt((0.25 * gamma_in) ** 2 - (Omega_in ** 2))
    g2_f = 1.0 - exp(-(3/4) * gamma_in * tau_in) * (cosh(d_f * tau_in) + \
                                                    ((3/4) * gamma_in / d_f) * sinh(d_f * tau_in))
    return real(g2_f)

def g2_dressed_states(tau_in, Omega_in, gamma_in, w0_in):
    """
    Analytic expression for the normalised second-order correlation function
    of the dressed state transitions.
    """
    from numpy import exp, ones
    if abs(w0_in) == abs(Omega_in):
        g2_f = 1.0 - exp(-0.5 * gamma_in * tau_in)
    elif abs(w0_in) == 0.0:
        g2_f = ones(len(tau_in))
    return g2_f

def pi_int(number_in):
    """
    Checks if the number is an integer multiple of pi and, if so, returns a 
    string r'{} \pi'.
    """
    from numpy import pi
    # If number_in is pretty much pi, set string to r'\pi'
    if round(number_in, 5) == round(pi, 5):
        str_out = r'\pi'

    # For larger-than-pi numbers
    if number_in > pi:
        if round(number_in % pi, 5) == 0:
            # Modulus is zero so number_in is an integer of pi
            # Set string
            str_out = r'{} \pi'.format(int(number_in // pi))
        else:
            # It's not an integer of pi so round it to 3dp
            str_out = r'{} \pi'.format(round(number_in, 3))
        
    # For less-than-pi numbers
    if number_in < pi:
        if round(pi % number_in, 5) == 0:
            # Modulus is zero so number_in is a rational of pi
            str_out = r'\frac{{\pi}}{{{}}}'.format(int(pi // number_in))
        else:
            # It's not an integer of pi so round it to 3dp
            str_out = r'{} \pi'.format(round(number_in, 3))
            
    # Return string
    return str_out

#-----------------------------------------------------------------------------#
#                               FILENAME THINGS                               #
#-----------------------------------------------------------------------------#
# Read parameters
gamma, Omega, N, halfwidth, kappa, dw, w0 = \
    np.genfromtxt("./data_files/g2_auto_parameters.txt", delimiter="=", usecols=1)
N = int(N)
phase = 1
# Omega = round(Omega, 2)
# w0 = round(w0, 2)

Omega_str = pi_int(Omega)
w0_str = pi_int(w0)

# Pull data from file
# tau
tau = np.genfromtxt("./data_files/g2_auto_corr.txt", usecols=0)
# g2
corr_filter = np.genfromtxt("./data_files/g2_auto_corr.txt", usecols=1)

# analytic
corr_atom = g2_analytic(tau, Omega, gamma)
# Dressed states g2
corr_dressed = g2_dressed_states(tau, Omega, gamma, w0)

print("Initial g2 value = {}".format(corr_filter[0]))

#-----------------------------------------------------------------------------#
#                               PLOT g^{(2)}                                  #
#-----------------------------------------------------------------------------#
fig = plt.figure(num='G2 Auto', figsize=[10, 6])
ax = plt.gca()

# Plot
if N > 0:
    ax.plot(tau, corr_filter, color='C0', ls='solid',
            label=(r'$N = {}, \delta\omega = {} \gamma, \kappa = {} \gamma$'
                   ).format(N, dw, kappa))
else:
    ax.plot(tau, corr_filter, color='C0', ls='solid',
            label='$\kappa = {} \gamma$'.format(kappa))

ax.plot(tau, corr_dressed, color='k', ls='dashed', alpha=0.5, label='Dressed-State Approximation')

# Legend
ax.legend(loc='lower right', fontsize=12)

# Labels
ax.set_xlabel(r'$\gamma \tau$', fontsize=12)
ax.set_ylabel(r'$g^{{(2)}}(\tau)$', fontsize=12)
ax.set_title((r'$g^{{(2)}}_{{\mathrm{{Auto}}}}(\tau)$ with $\Omega = {} \gamma'
              r', \omega_{{0}} = {} \gamma$'
              ).format(Omega_str, w0_str), fontsize=12)

# Figure stuff
fig.tight_layout()
fig.show()
