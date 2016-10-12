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
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    fi = FigInfo()

    # Plot settings
    st_fs = fi.t_fs              # Super Title Font Size
    sp_fs = fi.sp_fs             # Subplot Title Font Size

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
        save_name = fi.save_path + '201-CorrData' + '.' + fi.format
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


def plot_con_mat(dat, section, roi_lr, save_out=False):
    """Plot connectivity matrix.

    Parameters
    ----------
    dat : ?
        xx
    section : str
        Which part of the matrix to plot.
            Options: {'all', 'left', 'right', 'lr', 'rl'}
    roi_lr : list(str)
        List of L/R of ROIs.
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    fi = FigInfo()

    # Plot settings
    st_fs = fi.t_fs              # Super Title Font Size

    # Check the number of ROIs
    n_rois = len(roi_lr)

    # Get indices for section to plot
    ind_st_a, ind_en_a, ind_st_b, ind_en_b = get_section(section, n_rois, roi_lr)

    # Initialize figure
    f = plt.figure(figsize=(12, 12))
    ax = f.add_subplot(111)

    # Create the figure
    m = ax.imshow(dat[ind_st_a:ind_en_a, ind_st_b:ind_en_b], 
                  interpolation='none')

    # Add title
    plt.title('Connectivity Matrix', fontsize=st_fs, fontweight='bold')

    # Add colour bar as a legend
    f.colorbar(m)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + '202-ConnectivityMatrix' + '.' + fi.format
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)

