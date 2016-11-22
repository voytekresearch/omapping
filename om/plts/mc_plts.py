""" MODULE DOCSTRING - TO FILL IN
"""

from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np

from om.gen import FigInfo, get_section

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

