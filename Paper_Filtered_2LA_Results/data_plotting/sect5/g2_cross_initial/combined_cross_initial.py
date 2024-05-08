#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 14:27:57 2024

@author: jnga773
"""

import numpy as np
import matplotlib.pyplot as plt

# Add thesis style sheet
plt.style.use('../../paper_style.mplstyle')

plt.close('all')

# Filenames: Single-mode
par_single = "./data_files/SingleModeScan/parameters_initial.txt"
dat_single = "./data_files/SingleModeScan/g2_initial.txt"

# Filenames: Multi-mode
par_multi = "./data_files/MultiModeScan/parameters_initial.txt"
dat_multi = "./data_files/MultiModeScan/g2_initial.txt"

# Figure filename
filename_out = "../../../images/sect5/fig15_combined_cross_scan.pdf"

#-----------------------------------------------------------------------------#
#                                  FUNCTIONS                                  #
#-----------------------------------------------------------------------------#
def cmap_and_norm(boundaries_in, vcentre_in=1.0, ncolours_in=256,
                  top_colourmap='Blues_r', bottom_colourmap='Reds'):
    """
    Returns a diverging colourmap centred on vcentre_in and discrete boundary
    normalisation. The boundaries can be non-linear, and the colourmap
    will scale the non-linear boundaries automatically using 
    matplotlib.colors.BoundaryNorm.

    Instead of using one of matplotlib's default diverging colourmaps,
    the returned colourmap is made up of two separate colourmaps put on
    top of eachother, at the points vcentre_in.

    You can also set the centre-point of the colourmap, but it has to be
    one of the boundaries.
    
    The outputs, norm_out and cmap_out, can then be used in any of the
    matplotlib functions that need a colourmap and normalisation, e.g.
        plt.pcolormesh(Xdata, Ydata, Zdata, cmap=cmap_out, norm=norm_out)

    Parameters
    ----------
    boundaries_in : array-like
        A list (or array) of boundaries for the colourmap scale. Can be
        non-linear (e.g [0.01, 0.05, 0.2, 10.5, 125.23, 125.26, 3000.0]).
    vcentre_in : float (default = 1.0)
        The centre point of the diverging colourmap.
    ncolours_in : int (default = 256)
        Default number of RGB colours that is used in the colourmap. tbh
        I'm not entirely sure what this actually does haha
    top_colourmap : str (default = 'Blues_r')
        The colourmap colouring for the top half. Must be one of the
        colourmap names that matplotlib recognises.
    bottom_colourmap : str (default = 'Reds')
        The colourmap colouring for the bottom half. Must be one of the
        colourmap names that matplotlib recognises.
    Output
    ------
    norm_out : ???
        The normalisation used for the colourmap legend bar.
    colourmap_out : ???
        The colourmap for the data.
    """
    from numpy import interp, linspace, where, vstack
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from matplotlib.pyplot import get_cmap
    # First, check if vcentre is in the array boundaries_in
    for i in boundaries_in:
        if i == vcentre_in:
            vcentre_check = True
            break
        else:
            vcentre_check = False
    if vcentre_check is False:
        raise ValueError(
            "vcentre ({}) is not in boundaries_in!".format(vcentre_in))
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
    norm_cnst = norm_cnst[where(boundaries_in == vcentre_in)[0][0]]

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
    top = get_cmap(top_colourmap, no_top_colours)
    bottom = get_cmap(bottom_colourmap, no_bottom_colours)
    # top = colormaps.get_cmap(top_colourmap)
    # bottom = colormaps.get_cmap(bottom_colourmap)

    # Merge the colours maps
    new_colours = vstack((top(linspace(0, 1, no_top_colours)),
                         bottom(linspace(0, 1, no_bottom_colours))))
    cmap_out = ListedColormap(new_colours)

    return norm_out, cmap_out

#-----------------------------------------------------------------------------#
#                                  DATA THINGS                                #
#-----------------------------------------------------------------------------#
# Read parameters: Single-mode
gamma, Omega, N, K_single, kappa, dw  = \
    np.genfromtxt(par_single, delimiter="=", usecols=1, max_rows=6)

# Read parameters: Multi-mode
gamma, Omega, N, K_multi, kappa, dw  = \
    np.genfromtxt(par_multi, delimiter="=", usecols=1, max_rows=6)
N = int(N)

# Read data: Central freqency list
w_list = np.genfromtxt(par_multi, usecols=0, skip_header=8)

# Read data: Single-mode
g2_single = np.genfromtxt(dat_single, dtype='float')

# Read data: Multi-mode
g2_multi = np.genfromtxt(dat_multi, dtype='float')

# Create meshgrid of central frequency values
X, Y = np.meshgrid(w_list, w_list)

#-----------------------------------------------------------------------------#
#                           PLOT G2 SCAN: NON-LINEAR                          #
#-----------------------------------------------------------------------------#
# Boundary Norm
bounds = np.array([0.03, 0.1,
                   0.3, 1.0,
                   3.0, 10.0,
                   30.0, 100.0,
                   300.0, 1000.0,
                   3000.0])

# Create colourmap normalisation
norm, cmap = cmap_and_norm(bounds)

# Create string of bounds for labels
str_bounds = bounds
for place, i in enumerate(bounds):
    str_bounds[place] = str(i)

#--------------#
#     Plot     #
#--------------#
# from matplotlib import colorbar
fig, ax = plt.subplots(num='Initial Cross-Correlation Scan', nrows=1, ncols=2,
                        figsize=[12, 5])


# fig = plt.figure(num='Initial Cross-Correlation Scan', figsize=[12, 5],
#                  layout='compressed')
# from mpl_toolkits.axes_grid1 import ImageGrid
# ax = ImageGrid(fig, rect=111, nrows_ncols=(1, 2), axes_pad=0.15,
#                cbar_location='right', cbar_size="7%", cbar_pad=0.15)

# Plot: Single-mode
contour_single = ax[0].pcolormesh(X, Y, g2_single, shading='auto',
                                  cmap=cmap, norm=norm, rasterized=True)

# Plot: Multi-mode
contour_multi = ax[1].pcolormesh(X, Y, g2_multi, shading='auto',
                                 cmap=cmap, norm=norm, rasterized=True)

#--------------------#
#     Colour bar     #
#--------------------#
cbar_single = plt.colorbar(contour_single, ax=ax[0])
cbar_multi = plt.colorbar(contour_multi, ax=ax[1])

for cbar in [cbar_single, cbar_multi]:
    cbar.set_label(r'$g^{(2)}(\alpha, 0; \beta, 0)$', rotation=90)
    cbar.set_ticks(bounds)
    cbar.set_ticks([0.065, 0.2, 0.65, 2.0, 6.5, 20.0, 65.0, 200.0, 650.0, 2000.0], minor=True)
    
    cbar.ax.set_yticklabels(str_bounds)
    cbar.ax.tick_params(labelsize=8, direction='out', which='both')

#-------------------------#
#     Add Text Labels     #
#-------------------------#
# Label (a): Single-mode
ax[0].text(x=-45, y=32, s='(a)', fontsize=11)
# Label (b): Multi-mode
ax[1].text(x=-45, y=32, s='(b)', fontsize=11)

#--------------------#
#     Axis ticks     #
#--------------------#
w_ticks_maj = [-30, -Omega, 0.0, Omega, 30]
w_ticks_maj_str = [r'-$2\Omega$', r'$-\Omega$', r'0', r'$\Omega$', r'$2\Omega$']
w_ticks_min = [-1.5 * Omega, -0.5 * Omega, 0.5 * Omega, 1.5 * Omega]

for i in range(2):
    ax[i].set_xticks(w_ticks_maj, minor=False)
    ax[i].set_xticks(w_ticks_min, minor=True)
    
    ax[i].set_yticks(w_ticks_maj, minor=False)
    ax[i].set_yticks(w_ticks_min, minor=True)
    
    ax[i].set_xticklabels(w_ticks_maj_str, minor=False)
    ax[i].set_yticklabels(w_ticks_maj_str, minor=False)

#--------------#
#     Grid     #
#--------------#
for i in range(2):
    ax[i].tick_params(grid_alpha=0.75)
    ax[i].grid(color='lightgrey')

#--------------------#
#     Set labels     #
#--------------------#
for i in range(2):
    ax[i].set_xlabel(r'$\Delta\omega_{0}^{(a)}$')
    ax[i].set_ylabel(r'$\Delta\omega_{0}^{(b)}$')

#-----------------------#
#     Figures stuff     #
#-----------------------#
# Adjust subplot spacing
# plt.subplots_adjust(wspace=0.0)

# fig.tight_layout()
fig.savefig(filename_out)
fig.show()