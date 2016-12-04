""" MODULE DOCSTRING - TO FILL IN"""

from __future__ import print_function, division

import matplotlib.pyplot as plt
import numpy as np
import random

from om.plts.fig_info import FigInfo
from om.core.utils import get_section

##############################################################################################
############################## OM - PLTS - MAP COMPARISON PLOTS ##############################
##############################################################################################

def plot_corrs(corrs, p_vals, save_out=False):
    """Plot the R and p-values of the correlation results.

    Parameters
    ----------
    corrs : 1d array
        Vector of R values to plot distribution of.
    p_vals : 1d array
        Vector of R values to plot distribution of.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    f_info = FigInfo()

    # Plot settings
    st_fs = f_info.t_fs              # Super Title Font Size
    sp_fs = f_info.sp_fs             # Subplot Title Font Size

    # Check that input data is right dimension
    corrs = np.squeeze(corrs)
    p_vals = np.squeeze(p_vals)

    # Initialize subplot figure
    fig, ax = plt.subplots(1, 2, figsize=(12, 6))

    # Set super title for the plot
    plt.suptitle('Corr Stats', fontsize=st_fs, fontweight='bold')

    # Plot R values
    ax[0].plot(corrs, '.')
    ax[0].set_title('R Values', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[0].set_ylim([-0.5, 0.5])

    # Plot the p-values
    ax[1].plot(p_vals, '.')
    ax[1].set_title('P values', {'fontsize': sp_fs, 'fontweight': 'bold'})

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '201-CorrData' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


def plot_con_mat(dat, section, roi_lr, save_out=False):
    """Plot connectivity matrix.

    Parameters
    ----------
    dat : ?
        xx
    section : {'all', 'left', 'right', 'lr', 'rl'}
        Which part of the matrix to plot.
    roi_lr : list of str
        List of L/R of ROIs.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    f_info = FigInfo()

    # Plot settings
    st_fs = f_info.t_fs              # Super Title Font Size

    # Check the number of ROIs
    n_rois = len(roi_lr)

    # Get indices for section to plot
    ind_st_x, ind_en_x, ind_st_y, ind_en_y = get_section(section, n_rois, roi_lr)

    # Initialize figure
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111)

    # Create the figure
    mat = ax.imshow(dat[ind_st_x:ind_en_x, ind_st_y:ind_en_y],
                    interpolation='none')

    # Add title
    plt.title('Connectivity Matrix', fontsize=st_fs, fontweight='bold')

    # Add colour bar as a legend
    fig.colorbar(mat)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '202-ConnectivityMatrix' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


def plot_gene_corr(gene_dat, meg_dat, r_val, save_out=False):
    """Plot the gene & meg data for a particular gene/osc pair.

    Parameters
    ----------
    gene_dat : 1d array
        Vector of gene expression data.
    meg_dat : 1d array
        Vector of meg data to plot.
    r_val : float
        The actual
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    f_info = FigInfo()

    # Plot settings
    lw = f_info.ax_lw
    ax_fs = f_info.ax_fs
    alpha_val = 0.18

    # Extract only non-nan data
    inds_non_nan = np.invert(np.isnan(gene_dat))
    gene_dat = gene_dat[inds_non_nan]
    meg_dat = meg_dat[inds_non_nan]

    # Set number of points to plt
    n_points = 1000

    # Get a random sample of points to plot
    inds = random.sample(range(len(gene_dat)), n_points)

    # Initialize plot
    fig, ax = plt.subplots()

    # Draw scatter plot
    ax.scatter(gene_dat[inds], meg_dat[inds], color='#173570', marker='o', alpha=alpha_val)

    # Set axis limits
    space_meg = 0.05 * (max(meg_dat) - min(meg_dat))
    space_gene = 0.05 * (max(gene_dat) - min(gene_dat))

    # Figure out axis limits
    min_x = min(gene_dat) - space_gene
    max_x = max(gene_dat) + space_gene
    min_y = min(meg_dat) - space_meg
    max_y = max(meg_dat) + space_meg

    # Set axis limits on the plot
    plt.xlim([min_x, max_x])
    plt.ylim([min_y, max_y])

    # Add axis labels
    plt.xlabel('Gene Expression', {'fontsize': ax_fs})
    plt.ylabel('Oscillation Score', {'fontsize': ax_fs})

    # Set the top and right side frame & ticks off
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Set linewidth of remaining spines
    ax.spines['left'].set_linewidth(lw)
    ax.spines['bottom'].set_linewidth(lw)

    # Turn off other axis ticks
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    # Add text to the figure to report the r-value
    x_range = max_x - min_x
    ax.text((max_x-0.90*x_range), (3/4*max_y), 'r = ' + '{:4.4f}'.format(r_val),
            {'fontsize': 18, 'fontweight':'bold'})

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = cur_band + '_' + gene_name + '.pdf'
        plt.savefig(save_name, format='pdf', bbox_inches='tight', dpi=300)
