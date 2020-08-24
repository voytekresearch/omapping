"""OM - Plots for maps data."""

import random

import numpy as np
import matplotlib.pyplot as plt

from om.core.utils import get_section
from om.plts.fig_info import FigInfo
from om.plts.utils import set_axis_spines, save_figure

###################################################################################################
###################################################################################################

def plot_corrs(corrs, p_vals, save_out=False, fig_info=FigInfo()):
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

    # Plot settings
    st_fs = fig_info.t_fs              # Super Title Font Size
    sp_fs = fig_info.sp_fs             # Subplot Title Font Size

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

    save_figure(save_out, '201-CorrData', fig_info)


def plot_con_mat(data, section, roi_lr, save_out=False, fig_info=FigInfo()):
    """Plot connectivity matrix.

    Parameters
    ----------
    data : ?
        xx
    section : {'all', 'left', 'right', 'lr', 'rl'}
        Which part of the matrix to plot.
    roi_lr : list of str
        List of L/R of ROIs.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Plot settings
    st_fs = fig_info.t_fs              # Super Title Font Size

    # Check the number of ROIs
    n_rois = len(roi_lr)

    # Get indices for section to plot
    ind_st_x, ind_en_x, ind_st_y, ind_en_y = get_section(section, n_rois, roi_lr)

    # Initialize figure
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111)

    # Create the figure
    mat = ax.imshow(data[ind_st_x:ind_en_x, ind_st_y:ind_en_y],
                    interpolation='none')

    # Add title
    plt.title('Connectivity Matrix', fontsize=st_fs, fontweight='bold')

    # Add colour bar as a legend
    fig.colorbar(mat)

    save_figure(save_out, '202-ConnectivityMatrix', fig_info)


def plot_gene_corr(gene_data, meg_dat, r_val, save_out=False, fig_info=FigInfo()):
    """Plot the gene & meg data for a particular gene/osc pair.

    Parameters
    ----------
    gene_data : 1d array
        Vector of gene expression data.
    meg_data : 1d array
        Vector of meg data to plot.
    r_val : float
        The actual
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Plot settings
    lw = fig_info.ax_lw
    ax_fs = fig_info.ax_fs
    alpha_val = 0.18

    # Extract only non-nan data
    inds_non_nan = np.invert(np.isnan(gene_data))
    gene_data = gene_data[inds_non_nan]
    meg_data = meg_data[inds_non_nan]

    # Set number of points to plt
    n_points = 1000

    # Get a random sample of points to plot
    inds = random.sample(range(len(gene_data)), n_points)

    # Initialize plot
    fig, ax = plt.subplots()

    # Draw scatter plot
    ax.scatter(gene_data[inds], meg_data[inds], color='#173570', marker='o', alpha=alpha_val)

    # Set axis limits
    space_meg = 0.05 * (max(meg_data) - min(meg_data))
    space_gene = 0.05 * (max(gene_data) - min(gene_data))

    # Figure out axis limits
    min_x = min(gene_data) - space_gene
    max_x = max(gene_data) + space_gene
    min_y = min(meg_data) - space_meg
    max_y = max(meg_data) + space_meg

    # Set axis limits on the plot
    ax.set_xlim([min_x, max_x])
    ax.set_ylim([min_y, max_y])

    # Add axis labels
    ax.set_xlabel('Gene Expression', {'fontsize': ax_fs})
    ax.set_ylabel('Oscillation Score', {'fontsize': ax_fs})

    # Set spines
    set_axis_spines(ax, lw=fig_info.ax_lw)

    # Turn off other axis ticks
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])

    # Add text to the figure to report the r-value
    x_range = max_x - min_x
    ax.text((max_x-0.90*x_range), (3/4*max_y), 'r = ' + '{:4.4f}'.format(r_val),
            {'fontsize': 18, 'fontweight':'bold'})

    save_figure(save_out, cur_band + '_' + gene_name, fig_info)
