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

#-----------------------------------------------------------------------------#
#                                FILENAME THINGS                              #
#-----------------------------------------------------------------------------#
# Read parameters
gamma, Omega, kappa, dw, epsilon, N, phase = \
    np.genfromtxt(filename("parameters"), delimiter="=", skip_header=1, usecols=1)
N = int(N)
phase = int(phase)
Omega = round(Omega, 2)

# Pull data from file
# Central frequencies
w0_list = np.genfromtxt(filename("w0_cross"), usecols=0)
# Steady state photon number
g2_0 = np.genfromtxt(filename("cross_corr"))

X, Y =np.meshgrid(w0_list, w0_list)

#-----------------------------------------------------------------------------#
#                             PLOT G2 SCAN: LINEAR                            #
#-----------------------------------------------------------------------------#
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize, TwoSlopeNorm
# Set colormap
cmap = get_cmap('bwr', 256)

# Set normalisation
# norm = Normalize(vmin=g2_0.min(), vmax=g2_0.max())
norm = TwoSlopeNorm(vmin=g2_0.min(), vcenter=1.0, vmax=g2_0.max())

# Set up figure
fig, ax = plt.subplots(1, 1, figsize=[8, 7])

# Plot contourf
# contour_plot = ax.contourf(w0_list, w0_list, g2_0, 500,
#                             cmap=cmap, norm=norm)
contour_plot = ax.pcolormesh(X, Y, g2_0, shading='auto',
                              cmap=cmap, norm=norm)

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
# bounds = np.array([0.1, 0.3, 1.0, 3.0, 10.0, 30.0])
bounds = np.array([0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0])
norm, cmap = cmap_and_norm(bounds)

# Set up figure
fig, ax = plt.subplots(1, 1, figsize=[8, 7])

# Plot contourf
# contour_plot = ax.contourf(w0_list, w0_list, g2_0, 500,
#                             cmap=cmap, norm=norm)
contour_plot = ax.pcolormesh(X, Y, g2_0, shading='auto',
                              cmap=cmap, norm=norm)

# Colour bar
cbar = fig.colorbar(contour_plot, ax=ax, ticks=bounds,
                    label=r'$g^{(2)}_{\mathrm{Cross}}(\tau=0)$')

# Set labels and title
ax.set_xlabel(r'$\omega_{0}^{(a)} / \gamma$', fontsize=12)
ax.set_ylabel(r'$\omega_{0}^{(b)} / \gamma$', fontsize=12)
# ax.set_title(r'Initial Cross Correlation Value $g^{(2)}_{\mathrm{Cross}}(\tau=0)$', fontsize=12)
ax.set_title(r'Atom: $\Omega = {} \gamma$, Filter: $\kappa = {} \gamma, N = {}, \delta\omega = {} \gamma$'.format(Omega, kappa, N, dw), fontsize=12)

fig.tight_layout()
fig.show()