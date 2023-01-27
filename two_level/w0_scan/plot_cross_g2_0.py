# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:27:11 2019

@author: Jacob
"""

from _my_functions import filename

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'

plt.close('all')

#-----------------------------------------------------------------------------#
#                                FILENAME THINGS                              #
#-----------------------------------------------------------------------------#
# Read parameters
gamma, Omega, kappaa, dwa, kappab, dwb, epsilon, N, phase = \
    np.genfromtxt(filename("cross_parameters"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
Omega_str = pi_int(Omega)
w0_str = pi_int(w0)

# Pull data from file
# Central frequencies
w0a_list = np.genfromtxt(filename("w0_cross"), usecols=0)
w0b_list = np.genfromtxt(filename("w0_cross"), usecols=1)
# Steady state photon number
g2_0 = np.genfromtxt(filename("cross_corr"))

X, Y =np.meshgrid(w0a_list, w0b_list)

#-----------------------------------------------------------------------------#
#                             PLOT G2 SCAN: LINEAR                            #
#-----------------------------------------------------------------------------#
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize, TwoSlopeNorm

# Set colormap
if g2_0.min() < 1.0:
    cmap = get_cmap('bwr', 256)
    norm = TwoSlopeNorm(vmin=g2_0.min(), vcenter=1.0, vmax=g2_0.max())
else:
    cmap = get_cmap('Reds', 256)
    norm = Normalize(vmin=1.0, vmax=g2_0.max())

if g2_0.max() <= 1.0:
    cmap = get_cmap('Blues_r', 256)
    norm = Normalize(vmin=g2_0.min(), vmax=1.0)
    
# Print minimum and maximum
print("g^2 minimum: {}".format(g2_0.min()))
print("g^2 maximum: {}".format(g2_0.max()))

# Set normalisation

# Set up figure
fig, ax = plt.subplots(1, 1, figsize=[8, 7])

# Plot contourf
# contour_plot = ax.contourf(w0a_list, w0b_list, g2_0, 500,
#                             cmap=cmap, norm=norm)
contour_plot = ax.pcolormesh(X, Y, g2_0, shading='auto',
                              cmap=cmap, norm=norm)

# Plot cross
peaks = [-Omega, 0.0, Omega]
for i in peaks:
    for j in peaks:
        ax.plot(i, j, marker='x', color='green', lw=5.0)

# Colour bar
cbar = fig.colorbar(contour_plot, ax=ax,
                    label=r'$g^{(2)}_{\mathrm{Cross}}(\tau=0)$')

# Set labels and title
ax.set_xlabel(r'$\omega_{0}^{(a)} / \gamma$', fontsize=12)
ax.set_ylabel(r'$\omega_{0}^{(b)} / \gamma$', fontsize=12)
# ax.set_title(r'Initial Cross Correlation Value $g^{(2)}_{\mathrm{Cross}}(\tau=0)$', fontsize=12)

fig.tight_layout()
fig.show()

#-----------------------------------------------------------------------------#
#                           PLOT G2 SCAN: NON-LINEAR                          #
#-----------------------------------------------------------------------------#
#-----------------------------------------#
#     Set new colourmap/normalisation     #
#-----------------------------------------#
def cmap_and_norm(boundaries_in, vcentre=1.0, ncolours_in=256):
    """
    Returns a diverging colourmap centred on vcentre and discrete boundary
    normalisation.
    """
    from numpy import interp, linspace, where, vstack
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.cm import get_cmap
    # First, check if vcentre is in the array boundaries_in
    for i in boundaries_in:
        if i == vcentre:
            vcentre_check = True
            break
        else:
            vcentre_check = False
    if vcentre_check is False:
        raise ValueError("vcentre ({}) is not in boundaries_in!".format(vcentre))
    # Number of boundary points
    N_bounds = len(boundaries_in)
    # Stretch the bounds so colourmap is continuous-esque
    stretched_bounds = interp(linspace(0, 1, ncolours_in+1), linspace(0, 1, N_bounds),
                              boundaries_in)
    # Return the norm
    norm_out = BoundaryNorm(stretched_bounds, ncolors=ncolours_in)
    
    # The normalisations are taken from the boundaries_in
    norm_cnst = linspace(0, 1.0, N_bounds)
    # Find where vcentre is and take that fraction of colours to be the min
    norm_cnst = norm_cnst[where(boundaries_in == vcentre)[0][0]]
    
    # The colour map will constant of a reverse Blue on the bottom, merging
    # with a Red map on the top. The bottom will take up [norm_cnst]% of the
    # total colours [ncolours_in]
    no_top_colours = int(ncolours_in * norm_cnst)
    no_bottom_colours = int(ncolours_in * (1 - norm_cnst))
    # If the don't add up to ncolours_in, add some more until they do
    if no_top_colours + no_bottom_colours != ncolours_in:
        diff = ncolours_in - (no_top_colours + no_bottom_colours)
        no_top_colours += diff
    
    # Grab top and bottom colour maps
    top = get_cmap('Blues_r', no_top_colours)
    bottom = get_cmap('Reds', no_bottom_colours)
    
    # Merge the colours maps
    new_colours = vstack((top(linspace(0, 1, no_top_colours)), bottom(linspace(0, 1, no_bottom_colours))))
    cmap_out = ListedColormap(new_colours)
    
    return norm_out, cmap_out

#--------------#
#     Plot     #
#--------------#

# Boundary Norm
bounds = np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0])
# bounds = np.array([0.06, 0.1, 0.6, 1.0, 6.0, 10.0, 60.0, 100.0, 160.0])
# bounds = np.array([0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0])
norm, cmap = cmap_and_norm(bounds)

# Create string of bounds for labels
str_bounds = bounds
for place, i in enumerate(bounds):
    str_bounds[place] = str(i)

# Set up figure
fig, ax = plt.subplots(1, 1, figsize=[8, 7])

# Plot contourf
# contour_plot = ax.contourf(w0a_list, w0b_list, g2_0, 500,
#                             cmap=cmap, norm=norm)
contour_plot = ax.pcolormesh(X, Y, g2_0, shading='auto',
                              cmap=cmap, norm=norm)
# contour_plot = ax.imshow(g2_0, cmap=cmap, norm=norm)

# # Plot crosses on peaks
# peaks = [-Omega, 0.0, Omega]
# for i in peaks:
#     for j in peaks:
#         ax.plot(i, j, marker='x', color='k', lw=5.0)

# Colour bar
cbar = fig.colorbar(contour_plot, ax=ax, ticks=bounds,
                    label=r'$g^{(2)}_{\mathrm{Cross}}(\tau=0)$')
cbar.ax.set_yticklabels(str_bounds)
cbar.ax.tick_params(labelsize=15)

# Set ticks
w_ticks = [-30, -1.5 * Omega, -Omega, -0.5 * Omega, 0.0, 0.5 * Omega, Omega, 1.5 * Omega, 30]
w_ticks_str = [r'-$2\Omega$', r'$-1.5 \Omega$', r'$-\Omega$',
               r'$-0.5 \Omega$', '0', r'$0.5 \Omega$',
               r'$\Omega$', r'$1.5 \Omega$', r'$2\Omega$']
ax.set_xticks(w_ticks)
ax.set_yticks(w_ticks)
ax.set_xticklabels(w_ticks_str)
ax.set_yticklabels(w_ticks_str)

# Set labels and title
ax.set_xlabel(r'$\omega_{0}^{(a)} / \gamma$', fontsize=15)
ax.set_ylabel(r'$\omega_{0}^{(b)} / \gamma$', fontsize=15)
# ax.set_title(r'Initial Cross Correlation Value $g^{(2)}_{\mathrm{Cross}}(\tau=0)$', fontsize=12)
# ax.set_title(r'Atom: $\Omega = {} \gamma$, Filter: $\kappa = {} \gamma, N = {}, \delta\omega = {} \gamma$'.format(Omega, kappa, N, dw), fontsize=12)

fig.tight_layout()
fig.show()