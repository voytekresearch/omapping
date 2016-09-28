# Import required libraries/functions
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

#############################################################
#################### OM - PLTS - CLASSES ####################
#############################################################

class FigInfo():
    """Object to hold settings to save figures. """

    def __init__(self):

        # Default Settings
        self.t_fs = 20           # Title font size
        self.sp_fs = 18          # Subplot title font size
        self.ax_fs = 16          # Axis font size

        # Plot Information
        self.title = 'Group'
        self.vis_opac = 0.005

        # Save Information
        self.save_path = '/Users/thomasdonoghue/Desktop/'
        self.format = 'pdf'
        self.bbox = 'tight'
        self.dpi = 600

###############################################################
################## OM - PLTS - MEGDATA PLOTS ##################
###############################################################

def plot_slopes(slopes, title):
    """Plots a histogram of the chi values for all vertices. 

    Parameters
    ----------
    slopes : 1d array
        A vector of slope values to plot.
    title : str
        A string to append to the title.
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


def plot_hist_count(osc_count):
    """Plots a histogram of the osc_count vector.
    Parameters
    ----------
    osc_count : 1d Vector
        An array of the number of oscillations found in each vertex. 
    """

    # Get FigInto
    fi = FigInfo()

    # Plot Settings
    n_bins = 25              # Number of bins for histograms
    t_fs = fi.t_fs           # Title font size
    ax_fs = fi.ax_fs         # Axis label font size

    # Create histogram
    plt.hist(osc_count, n_bins, range=[0, 8])
    plt.title('# Oscillations per Vertex', {'fontsize': t_fs, 'fontweight': 'bold'})
    plt.xlabel('# Oscillations', {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})


def plot_all_oscs(centers_all, powers_all, bws_all, title):
    """Plots histogram distributions of oscillation centers, powers and bws. """

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


def plot_comparison(centers_all, powers_all, bws_all, title):
    """Plots comparisons between oscillatory parameters.

    Checks Centers vs. Bandwidth, Centers vs. Power and Bandwidth vs. Power.
    """

    # Get FigInto
    fi = FigInfo()

    # Plot Settings
    st_fs = fi.t_fs              # Super Title Font Size
    sp_fs = fi.sp_fs             # Subplit Title Font Size
    ax_fs = fi.ax_fs             # Axis Label Font Size
    vis_opac = 0.1              # Alpha value for plotted data

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


####################################################################
################## OM - PLTS - GROUPMEGDATA PLOTS ##################
####################################################################


def plot_freq_corr(fs, corr_vec, p_vec):
    """Creats a scatter plot for the rolling frequency correlation.

    Parameters
    ----------
    fs : 1d array
        Vector of frequencies used for rolling correlation (x-data).
    corr_vec : 1d array
        Vector of correlations between frequency bins (y-data).
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
    fig, ax = plt.subplots(figsize=(8,5))

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


def plot_age_peak(age, peak_theta, peak_alpha, peak_beta, peak_lowgamma):
    """   """

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


def plot_age_n_oscs(ages, nOscs):
    """Create a scatter plot comparing age and number of oscillations."""

    # Get FigInfo
    fi = FigInfo()

    # Plot settings
    t_fs = fi.t_fs
    ax_fs = fi.ax_fs

    # Make the plot
    plt.plot(ages, nOscs, '.')
    plt.title('# Oscillations / Age', fontsize=20, fontweight='bold')

