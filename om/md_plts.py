# Import required libraries/functions
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from om.gen import FigInfo

###############################################################################
########################## OM - PLTS - MEGDATA PLOTS ##########################
###############################################################################

def plot_slopes(slopes, title, save_out=False):
    """Plots a histogram of the chi values for all vertices.

    Parameters
    ----------
    slopes : 1d array
        A vector of slope values to plot.
    title : str
        A string to append to the title.
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInfo()
    fi = FigInfo()

    # Plot Settings
    n_bins = 100             # Number of bins for histograms
    t_fs = fi.t_fs           # Title font size
    ax_fs = fi.ax_fs         # Axis label font size

    # Create histogram
    plt.hist(slopes, n_bins)
    plt.title('Slopes - ' + title, {'fontsize': t_fs, 'fontweight': 'bold'})
    plt.xlabel('Chi Parameter', {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + title + '_Slopes.pdf'
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


def plot_hist_count(osc_count, save_out=False):
    """Plots a histogram of the osc_count vector.

    Parameters
    ----------
    osc_count : 1d Vector
        An array of the number of oscillations found in each vertex.
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    fi = FigInfo()

    # Plot Settings
    n_bins = 25              # Number of bins for histograms
    t_fs = fi.t_fs           # Title font size
    ax_fs = fi.ax_fs         # Axis label font size

    # Create histogram
    plt.hist(osc_count, n_bins, range=[0, 8])

    # Add title and axis labels
    plt.title('# Oscillations per Vertex', {'fontsize': t_fs, 'fontweight': 'bold'})
    plt.xlabel('# Oscillations', {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + 'OscCount.pdf'
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


def plot_all_oscs(centers_all, powers_all, bws_all, title, save_out=False):
    """Plots combined plot with distributions of oscillation centers, powers and bws.

    Parameters
    ----------
    centers_all : 1d array
        Vector of the center frequencies for all oscillations.
    powers_all : 1d array
        Vector of the powers for all oscillations.
    bws_all : 1d array
        Vector of the bandwidths for all oscillations.
    title : str
        A string to append to the title.
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    fi = FigInfo()

    # Plot Settings
    n_bins = 160             # Number of bins for histograms
    st_fs = fi.t_fs          # Super Title Font Size
    sp_fs = fi.sp_fs         # Subplot Title Font Size
    ax_fs = fi.ax_fs         # Axis Label Font Size

    # Set up subplots
    fig, ax = plt.subplots(3, 1, figsize=(15, 15))

    # Set plot super-title
    plt.suptitle('Distributions of Oscillatory Parameters - ' + title,
                 fontsize=st_fs, fontweight='bold')

    # Subplot 1 - Center Frequency
    ax[0].hist(centers_all, n_bins)
    ax[0].set_title('Center Frequency', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[0].set_xlabel('Frequency', {'fontsize': ax_fs})
    ax[0].set_ylabel('Count', {'fontsize': ax_fs})

    # Subplot 2 - Power
    ax[1].hist(np.log10(powers_all), n_bins)
    ax[1].set_title('Oscillatory Power', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[1].set_xlabel('Log Power', {'fontsize': ax_fs})
    ax[1].set_ylabel('Count', {'fontsize': ax_fs})

    # Subplot 3 - Bandwidth
    ax[2].hist(bws_all, n_bins)
    ax[2].set_title('Band Width', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[2].set_xlabel('Bandwidth (Hz)', {'fontsize': ax_fs})
    ax[2].set_ylabel('Count', {'fontsize': ax_fs})

    # Adjust subplot spacing
    plt.subplots_adjust(hspace=0.4)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + 'AllOscs.pdf'
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


def plot_all_osc_single(data, dat_type, title, n_bins=160, size=(5, 15), save_out=False):
    """Create a plot for a single oscillation parameter.

    Parameters
    ----------
    data : 1d array
        Vector of oscillation parameter data.
    dat_type : str
        Which osc parameter is being plotted.
    title : str
        A string to append to the title.
    n_bins : int, optional
        Number of bins to use for the plot. Default: 160.
    size : tuple, optional
        Size of the figure to make. Default: (5, 15).
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    fi = FigInfo()

    # Plot Settings
    t_fs = fi.t_fs           # Super Title Font Size
    ax_fs = fi.ax_fs         # Axis Label Font Size

    # Set up plot
    fig, ax = plt.subplots(figsize=size)

    # Subplot 1 - Center Frequency
    ax.hist(centers_all, n_bins)

    plt.title('Center Frequency', {'fontsize': sp_fs, 'fontweight': 'bold'})
    plt.xlabel('Frequency', {'fontsize': ax_fs})
    plt.ylabel('Count', {'fontsize': ax_fs})

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + dat_type + '_AllOscs.pdf'
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


def plot_osc_param_comparison(centers_all, powers_all, bws_all, title, save_out=False):
    """Plots comparisons between all oscillatory parameters.

    Checks Centers vs. Bandwidth, Centers vs. Power and Bandwidth vs. Power.

    Parameters
    ----------
    centers_all : 1d array
        Vector of the center frequencies for all oscillations.
    powers_all : 1d array
        Vector of the powers for all oscillations.
    bws_all : 1d array
        Vector of the bandwidths for all oscillations.
    title : str
        A string to append to the title.
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    fi = FigInfo()

    # Plot Settings
    st_fs = fi.t_fs              # Super Title Font Size
    sp_fs = fi.sp_fs             # Subplit Title Font Size
    ax_fs = fi.ax_fs             # Axis Label Font Size
    vis_opac = 0.1               # Alpha value for plotted data

    # Set up subplots
    fig, ax = plt.subplots(3, 1, figsize=(15, 15))

    # Set plot super-title
    plt.suptitle('Oscillation Parameter Comparisons - ' + title,
                 fontsize=st_fs, fontweight='bold')

    # Plot - Center vs. Bandwidth
    ax[0].plot(centers_all, np.log10(bws_all), '.', alpha=vis_opac)
    ax[0].set_title('Center vs. Bandwidth', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[0].set_xlabel('Centers', {'fontsize': ax_fs})
    ax[0].set_ylabel('BW', {'fontsize': ax_fs})

    # Plot - Centers vs. Power
    ax[1].plot(centers_all, np.log10(powers_all), '.', alpha=vis_opac)
    ax[1].set_title('Center vs. Power', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[1].set_xlabel('Centers', {'fontsize': ax_fs})
    ax[1].set_ylabel('Log Power', {'fontsize': ax_fs})

    # Plot - BWs vs. Powers
    ax[2].plot(np.log10(bws_all), np.log10(powers_all), '.', alpha=vis_opac)
    ax[2].set_title('BW vs. Power', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[2].set_xlabel('Bandwidth (Hz)', {'fontsize': ax_fs})
    ax[2].set_ylabel('Log Power', {'fontsize': ax_fs})

    # Adjust subplot spacing
    plt.subplots_adjust(hspace=0.4)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + 'OscComparison.pdf'
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


########################################################################################
############################ OM - PLTS - GROUPMEGDATA PLOTS ############################
########################################################################################


def plot_peak_boxplot(peaks, osc_band, save_out=False):
    """Plot a boxplot of peak frequencies within an oscillation band.

    Parameters
    ----------
    peaks : 1d array
        Vector of peak frequencies for given oscillation band.
    osc_band : str
        Label of which osc band is being plotted.
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Initialize the figure
    f, ax = plt.subplots(figsize=(2, 5))

    # Make the boxplot
    ax.boxplot(peaks, widths=0.45)

    # Set the top and right side frame & ticks off
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Set linewidth of remaining spines
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

    # Add a title
    plt.title(osc_band, fontsize=20, fontweight='bold')

    # Remove x ticks
    plt.xticks([])

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + osc_band + '_boxplot.pdf'
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


def plot_freq_corr(fs, corr_vec, p_vec, save_out=False):
    """Creats a scatter plot for the rolling frequency correlation.

    Parameters
    ----------
    fs : 1d array
        Vector of frequencies used for rolling correlation (x-data).
    corr_vec : 1d array
        Vector of correlations between frequency bins (y-data).
    p_vec : 1d array
        Vector of p-values for each correlation.
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInfo
    fi = FigInfo()

    # Plot settings
    t_fs = fi.t_fs
    ax_fs = fi.ax_fs
    r_alpha = 1
    r_size = 10
    p_alpha = 0.6
    p_size = 5

    # Make the plot
    fig, ax = plt.subplots(figsize=(8, 5))

    # Create the plot
    ax.plot(fs, corr_vec, 'x', alpha=r_alpha, markersize=r_size)

    # Add marker for significance
    for i, p in enumerate(p_vec):
        
        # Add a marker if passes Bonferroni corrected p-value
        if p < 0.05/len(p_vec):
            ax.plot(fs[i], corr_vec[i], 'ro', alpha=p_alpha, markersize=p_size)

    # Set title
    plt.suptitle('Correlation Adjacent Frequency Bands', fontsize=t_fs, fontweight='bold')

    # Set axis limits
    ax.set_xlim([0, 38])
    ax.set_ylim([-1.0, 1.0])

    # Set axis labels
    ax.set_xlabel('Frequency', fontsize=ax_fs, fontweight='bold')
    ax.set_ylabel('R Value', fontsize=ax_fs, fontweight='bold')

    # Set the top and right side frame & ticks off
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Set linewidth of remaining spines
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + 'FreqCorr.pdf'
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


def plot_age_peak(age, peak_theta, peak_alpha, peak_beta, peak_lowgamma, save_out=False):
    """Createa a plot comparing age to peak frequencies for each oscillation band.

    Parameters
    ----------
    age : 1d array
        Vector of ages for each subject.
    peak_theta : 1d array
        Vector of peak theta frequency for each subject.
    peak_alpha : 1d array
        Vector of peak alpha frequency for each subject.
    peak_beta : 1d array
        Vector of peak beta frequency for each subject.
    peak_lowgamma : 1d array
        Vector of peak low gamma frequency for each subject.
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInfo
    fi = FigInfo()

    # Plot settings
    st_fs = fi.t_fs
    sp_fs = fi.sp_fs

    # Set up subplots
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))

    # Set plot super-title
    plt.suptitle('Peak Frequency / Age Comparisons', fontsize=st_fs, fontweight='bold')

    # Theta
    ax[0, 0].plot(age, peak_theta, '.')
    ax[0, 0].set_title('Theta', {'fontsize': sp_fs, 'fontweight': 'bold'})

    # Alpha
    ax[0, 1].plot(age, peak_alpha, '.')
    ax[0, 1].set_title('Alpha', {'fontsize': sp_fs, 'fontweight': 'bold'})

    # Beta
    ax[1, 0].plot(age, peak_beta, '.')
    ax[1, 0].set_title('Beta', {'fontsize': sp_fs, 'fontweight': 'bold'})

    # Gamma
    ax[1, 1].plot(age, peak_lowgamma, '.')
    ax[1, 1].set_title('Low Gamma', {'fontsize': sp_fs, 'fontweight': 'bold'})

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + 'AgePeak.pdf'
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


def plot_age_n_oscs(ages, n_oscs, save_out=False):
    """Create a scatter plot comparing age and number of oscillations.

    Parameters
    ----------
    ages : 1d array
        Vector of ages for each subject.
    n_oscs : 1d array
        Vector of the number of oscillations for each subject.
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInfo
    fi = FigInfo()

    # Plot settings
    t_fs = fi.t_fs
    ax_fs = fi.ax_fs

    # Make the plot
    plt.plot(ages, n_oscs, '.')
    plt.title('# Oscillations / Age', fontsize=20, fontweight='bold')

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + 'AgeNumberOscillations.pdf'
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


def plot_osc_profiles(centers_hist, save_out=False):
    """Creates a plot showing all subjects oscillation profiles.

    Parameters
    ----------
    centers_hist : list[list]
        Contains the oscillation profile for each subject.
    save_out : boolean, optional
        Whether to save out a copy of the figure.
    """

    # Get FigInfo
    fi = FigInfo()

    # Plot settings
    t_fs = fi.t_fs
    ax_fs = fi.ax_fs
    ind_lw = 0.5
    avg_lw = 2
    alpha = 0.4

    # Initialize plot
    fig, ax = plt.subplots(figsize=(18, 5))

    # Initialize a frequency vector for the x-axis
    freqs = np.arange(3.125, 40.125, 0.25)
    
    # Loop through all subjects, adding profile to plot
    for h in centers_hist:
        ax.plot(freqs, h, linewidth=ind_lw, alpha=alpha)
    
    # Add the average profile
    ax.plot(freqs, np.median(centers_hist, 0), 'k', linewidth=avg_lw)

    # Add title and axis labels
    plt.title('Individual Oscillation Profiles', {'fontsize': t_fs, 'fontweight': 'bold'})
    plt.xlabel('Frequency', {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Set the x axis limits
    plt.xlim([3, 40])

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = fi.save_path + 'OscillationProfiles.pdf'
        plt.savefig(save_name, format=fi.format, bbox_inches=fi.bbox, dpi=fi.dpi)


#######################################################################################
############################ OM - PLTS - MAP COMPARE PLOTS ############################
#######################################################################################


