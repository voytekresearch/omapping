"""MODULE DOCSTING - TO FILL IN"""

# Import required libraries/functions
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

#from om.gen import FigInfo, UnknownDataTypeError
from om.plts.fig_info import FigInfo
from om.core.errors import UnknownDataTypeError

#######################################################################################
############################## OM - PLTS - MEGDATA PLOTS ##############################
#######################################################################################

def plot_slopes(slopes, title, save_out=False):
    """Plots a histogram of the chi values for all vertices.

    Parameters
    ----------
    slopes : 1d array
        A vector of slope values to plot.
    title : str
        A string to append to the title.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInfo()
    f_info = FigInfo()

    # Plot Settings
    n_bins = 100             # Number of bins for histograms
    t_fs = f_info.t_fs           # Title font size
    ax_fs = f_info.ax_fs         # Axis label font size
    ti_fs = f_info.ti_fs         # Axis ticks font size

    # Create histogram
    plt.hist(slopes, n_bins)

    # Add title
    if f_info.add_title:
        plt.title('Slopes - ' + title, {'fontsize': t_fs, 'fontweight': 'bold'})

    # Add axis labels
    plt.xlabel('Chi Parameter', {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '101-' + title + '_Slopes' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


def plot_hist_count(osc_count, save_out=False):
    """Plots a histogram of the osc_count vector.

    Parameters
    ----------
    osc_count : 1d Vector
        An array of the number of oscillations found in each vertex.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    f_info = FigInfo()

    # Plot Settings
    n_bins = 25              # Number of bins for histograms
    t_fs = f_info.t_fs           # Title font size
    ax_fs = f_info.ax_fs         # Axis label font size
    ti_fs = f_info.ti_fs         # Axis ticks font size

    # Create histogram
    plt.hist(osc_count, n_bins, range=[0, 8])

    # Add title
    if f_info.add_title:
        plt.title('# Oscillations per Vertex', {'fontsize': t_fs, 'fontweight': 'bold'})

    # Add axis labels
    plt.xlabel('# Oscillations', {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '102-OscCount' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


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
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    f_info = FigInfo()

    # Plot Settings
    n_bins = 160             # Number of bins for histograms
    st_fs = f_info.t_fs          # Super Title Font Size
    sp_fs = f_info.sp_fs         # Subplot Title Font Size
    ax_fs = f_info.ax_fs         # Axis Label Font Size
    ti_fs = f_info.ti_fs         # Axis ticks font size

    # Set up subplots
    fig, ax = plt.subplots(3, 1, figsize=(15, 15))

    # Set plot super-title
    if f_info.add_title:
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

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Adjust subplot spacing
    plt.subplots_adjust(hspace=0.4)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '103-AllOscs' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


def plot_all_oscs_single(data, dat_type, title, n_bins=160, size=(15, 5), save_out=False):
    """Create a plot for a single oscillation parameter.

    Parameters
    ----------
    data : 1d array
        Vector of oscillation parameter data.
    dat_type : {0, 1, 2}
        Int refers to which osc parameter is being plotted.
            Key: {0:'Center Frequency', 1:'Power', 2:'Bandwidth'
    title : str
        A string to append to the title.
    n_bins : int, optional (default = 160)
        Number of bins to use for the plot.
    size : tuple, optional (default = (15, 5))
        Size of the figure to make.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    f_info = FigInfo()

    # Plot Settings
    t_fs = f_info.t_fs           # Super Title Font Size
    ax_fs = f_info.ax_fs         # Axis Label Font Size
    ti_fs = f_info.ti_fs         # Axis ticks font size
    ax_lw = f_info.ax_lw

    # Set up for which data type
    if dat_type is 0:
        dat_title = 'Center Frequency'
        xlab = 'Center Frequency'
    elif dat_type is 1:
        dat_title = 'Power'
        xlab = 'Log Power'
        data = np.log10(data)
    elif dat_type is 2:
        dat_title = 'Bandwidth'
        xlab = 'Bandwidth (Hz)'
    else:
        raise UnknownDataTypeError('Data type not understood.')

    # Set up plot
    fig, ax = plt.subplots(figsize=size)

    # Subplot 1 - Center Frequency
    ax.hist(data, n_bins, color='#40425e', alpha=0.9)

    # Set the title
    if f_info.add_title:
        plt.title(dat_title + ' ' + title, {'fontsize': t_fs, 'fontweight': 'bold'})

    # Add axis labels
    plt.xlabel(xlab, {'fontsize': ax_fs})
    plt.ylabel('Count', {'fontsize': ax_fs})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Set the top and right side frame & ticks off
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Set linewidth of remaining spines
    ax.spines['left'].set_linewidth(ax_lw)
    ax.spines['bottom'].set_linewidth(ax_lw)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '104-' + dat_title + '_AllOscs'+ '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


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
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    f_info = FigInfo()

    # Plot Settings
    st_fs = f_info.t_fs              # Super Title Font Size
    sp_fs = f_info.sp_fs             # Subplot Title Font Size
    ax_fs = f_info.ax_fs             # Axis Label Font Size
    ti_fs = f_info.ti_fs             # Axis ticks font size
    vis_opac = 0.1                   # Alpha value for plotted data

    # Set up subplots
    fig, ax = plt.subplots(3, 1, figsize=(15, 15))

    # Set plot super-title
    if f_info.add_title:
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

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '105-OscComparison' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


########################################################################################
############################ OM - PLTS - GROUPMEGDATA PLOTS ############################
########################################################################################


def plot_band_corr_matrix(corr_dat, save_out=False):
    """Plot the correlation ...

    Parameters
    ----------
    corr_dat : ?
        xx
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    f_info = FigInfo()

    # Plot Settings
    t_fs = f_info.t_fs              # Title Font Size
    ti_fs = f_info.ti_fs            # Axis ticks font size

    # Set colormap to use
    cmap = plt.get_cmap('seismic')

    # Create the plot
    #im = plt.matshow(corr_dat, vmin=-1, vmax=1, cmap=cmap, interpolation='none')
    im = plt.matshow(corr_dat, vmin=-1, vmax=1, cmap=cmap, interpolation='nearest')
    # Notes on using nearest here:
    #   https://github.com/matplotlib/matplotlib/issues/2972/

    # Add title
    if f_info.add_title:
        plt.title('Osc Band Correlations', {'fontsize': t_fs, 'fontweight': 'bold'}, y=1.15)

    # Set tick labels
    plt.xticks([0, 1, 2, 3], ['Theta', 'Alpha', 'Beta', 'LowGamma'], rotation=45, ha='left')
    plt.yticks([0, 1, 2, 3], ['Theta', 'Alpha', 'Beta', 'LowGamma'])

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Add a colorbar - add padding to offset further from plot
    plt.colorbar(pad=0.15)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '106-CorrMat' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


def plot_peak_boxplot(peaks, osc_band, save_out=False):
    """Plot a boxplot of peak frequencies within an oscillation band.

    Parameters
    ----------
    peaks : 1d array
        Vector of peak frequencies for given oscillation band.
    osc_band : str
        Label of which osc band is being plotted.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInto
    f_info = FigInfo()

    # Plot Settings
    t_fs = f_info.sp_fs              # Title Font Size
    ti_fs = f_info.ti_fs            # Axis ticks font size
    ax_lw = f_info.ax_lw

    # Initialize the figure
    fig, ax = plt.subplots(figsize=(2, 5))

    # Make the boxplot
    ax.boxplot(peaks, widths=0.45)

    # Set the top and right side frame & ticks off
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Set linewidth of remaining spines
    ax.spines['left'].set_linewidth(ax_lw)
    ax.spines['bottom'].set_linewidth(ax_lw)

    # Add a title
    #if fi.add_title:
    plt.title(osc_band, fontsize=t_fs, fontweight='bold', y=1.08)

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Remove x ticks
    plt.xticks([])

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '107-' + osc_band + '_boxplot' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


def plot_peak_boxplot_all(peaks, bands, save_out=False):
    """Plot boxplots of peak frequencies for all oscillation bands.

    Parameters
    ----------
    peaks : ?
        xx
    bands : ?
        xx

    TODO: Figure out how to space out subplots a little
    """

    # Get FigInto
    f_info = FigInfo()

    # Plot Settings
    t_fs = f_info.sp_fs             # Title Font Size
    ti_fs = f_info.ti_fs            # Axis ticks font size
    ax_lw = f_info.ax_lw            # Axis line weight

    n_bands = len(bands)
    f, ax = plt.subplots(1, n_bands, figsize=(10, 6))

    for ind, band in enumerate(bands):
        ax[ind].boxplot(peaks[band], widths=0.40)

        # Add title and set y lims using bands and data
        ax[ind].set_title(band, fontsize=t_fs, fontweight='bold', y=1.04)
        ax[ind].set_ylim([peaks[band].min()-0.2, peaks[band].max()+0.2])

    for ai in ax:

        # Set the top and right side frame & ticks off
        ai.spines['right'].set_visible(False)
        ai.spines['top'].set_visible(False)
        ai.xaxis.set_ticks_position('bottom')
        ai.yaxis.set_ticks_position('left')

        # Set linewidth of remaining spines
        ai.spines['left'].set_linewidth(ax_lw)
        ai.spines['bottom'].set_linewidth(ax_lw)

        # Set ticks font size
        ai.tick_params(axis='both', which='major', labelsize=ti_fs)

        # Remove x ticks
        ai.set_xticks([])

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + 'XXX-PeakFreqs_All' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


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
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInfo
    f_info = FigInfo()

    # Plot settings
    t_fs = f_info.t_fs
    ax_fs = f_info.ax_fs
    ti_fs = f_info.ti_fs         # Axis ticks font size
    ax_lw = f_info.ax_lw

    #
    r_alpha = 1
    r_size = 10
    p_alpha = 0.6
    p_size = 5

    # Make the plot
    fig, ax = plt.subplots(figsize=(8, 5))

    # Create the plot
    ax.plot(fs, corr_vec, 'x', alpha=r_alpha, markersize=r_size)

    # Add marker for significance
    for i, p_val in enumerate(p_vec):

        # Add a marker if passes Bonferroni corrected p-value
        if p_val < 0.05/len(p_vec):
            ax.plot(fs[i], corr_vec[i], 'ro', alpha=p_alpha, markersize=p_size)

    # Set title
    if f_info.add_title:
        plt.suptitle('Correlation Adjacent Frequency Bands', fontsize=t_fs, fontweight='bold')

    # Set axis limits
    ax.set_xlim([0, 38])
    ax.set_ylim([-1.0, 1.0])

    # Set axis labels
    ax.set_xlabel('Frequency', fontsize=ax_fs, fontweight='bold')
    ax.set_ylabel('R Value', fontsize=ax_fs, fontweight='bold')

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Set the top and right side frame & ticks off
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Set linewidth of remaining spines
    ax.spines['left'].set_linewidth(ax_lw)
    ax.spines['bottom'].set_linewidth(ax_lw)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '108-FreqCorr' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


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
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInfo
    f_info = FigInfo()

    # Plot settings
    st_fs = f_info.t_fs
    sp_fs = f_info.sp_fs
    ti_fs = f_info.ti_fs         # Axis ticks font size

    # Set up subplots
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))

    # Set plot super-title
    if f_info.add_title:
        plt.suptitle('Peak Frequency / Age Comparisons', fontsize=st_fs, fontweight='bold')

    # Theta
    ax[0, 0].plot(age, peak_theta, '.')
    ax[0, 0].set_title('Theta', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[0, 0].tick_params(axis='both', which='major', labelsize=ti_fs)

    # Alpha
    ax[0, 1].plot(age, peak_alpha, '.')
    ax[0, 1].set_title('Alpha', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[0, 1].tick_params(axis='both', which='major', labelsize=ti_fs)

    # Beta
    ax[1, 0].plot(age, peak_beta, '.')
    ax[1, 0].set_title('Beta', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[1, 0].tick_params(axis='both', which='major', labelsize=ti_fs)

    # Gamma
    ax[1, 1].plot(age, peak_lowgamma, '.')
    ax[1, 1].set_title('Low Gamma', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[1, 1].tick_params(axis='both', which='major', labelsize=ti_fs)

    # Set ticks font size
    #plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '109-AgePeak' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


def plot_age_n_oscs(ages, n_oscs, save_out=False):
    """Create a scatter plot comparing age and number of oscillations.

    Parameters
    ----------
    ages : 1d array
        Vector of ages for each subject.
    n_oscs : 1d array
        Vector of the number of oscillations for each subject.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInfo
    f_info = FigInfo()

    # Plot settings
    t_fs = f_info.t_fs
    ti_fs = f_info.ti_fs         # Axis ticks font size
    ax_fs = f_info.ax_fs

    # Make the plot
    plt.plot(ages, n_oscs, '.')

    # Add title
    if f_info.add_title:
        plt.title('# Oscillations / Age', fontsize=t_fs, fontweight='bold')

    # Add axis labels
    plt.xlabel('Age', {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '110-AgeNumberOscillations' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


def plot_osc_profiles(centers_hist, save_out=False):
    """Creates a plot showing all subjects oscillation profiles.

    Parameters
    ----------
    centers_hist : list of list
        Contains the oscillation profile for each subject.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Get FigInfo
    f_info = FigInfo()

    # Plot settings
    t_fs = f_info.t_fs
    ax_fs = f_info.ax_fs
    ti_fs = f_info.ti_fs         # Axis ticks font size
    ax_lw = f_info.ax_lw

    # xx
    ind_lw = 0.5
    avg_lw = 2.5
    alpha = 0.45

    # Initialize plot
    fig, ax = plt.subplots(figsize=(18, 5))

    # Initialize a frequency vector for the x-axis
    freqs = np.arange(3.125, 40.125, 0.25)

    # Loop through all subjects, adding profile to plot
    for hist in centers_hist:
        ax.plot(freqs, hist, linewidth=ind_lw, alpha=alpha)

    # Add the average profile
    ax.plot(freqs, np.median(centers_hist, 0), 'k', linewidth=avg_lw)

    # Add title
    if f_info.add_title:
        plt.title('Individual Oscillation Profiles', {'fontsize': t_fs, 'fontweight': 'bold'})

    # Set the top and right side frame & ticks off
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Set linewidth of remaining spines
    ax.spines['left'].set_linewidth(ax_lw)
    ax.spines['bottom'].set_linewidth(ax_lw)

    # Add axis labels
    plt.xlabel('Frequency', {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Set the x axis limits
    plt.xlim([3, 40])

    # Save out (if requested)
    if save_out:

        # Set up save name & save out
        save_name = f_info.save_path + '111-OscillationProf_infoles' + '.' + f_info.format
        plt.savefig(save_name, format=f_info.format, bbox_inches=f_info.bbox, dpi=f_info.dpi)


def plot_space_comp(oscs, verts, band, subj=0, osc_param=0, space_param=1):
    """

    Parameters
    ----------
    oscs : dict()
        xx
    verts : ?
        xx
    band : ?
        xx
    subj : ?
        xx
    osc_param : ?
        xx
    space_param : ?
        xx
    """

    # Get FigInfo
    f_info = FigInfo()

    # Plot settings
    t_fs = f_info.t_fs
    ax_fs = f_info.ax_fs

    #
    space = verts[:, space_param]
    sort_inds = np.argsort(space)
    freqs = [dat if dat > 0 else None for dat in oscs[band][sort_inds, osc_param, subj]]

    # Initialize Figure
    fig = plt.figure()

    #
    plt.plot(space, freqs, '.', ms=3.5, alpha=0.75)

    # Set plotting limits
    plt.xlim([space.min()-4, space.max()+4])

    #
    plt.title(band, {'fontsize': t_fs, 'fontweight': 'bold'})

    #
    plt.xlabel('Posterior -> Anterior', {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel('Center Frequency', {'fontsize': ax_fs, 'fontweight': 'bold'})


def plot_space_comp_all():
    """   """
    pass


def plot_osc_space_corr_boxplot(dat, labels):
    """

    Parameters
    ----------
    dat : ?
        xx
    labels : ?
        xx
    """

    # Get FigInfo
    f_info = FigInfo()

    # Plot settings
    t_fs = f_info.t_fs
    ax_fs = f_info.ax_fs
    ti_fs = f_info.ti_fs         # Axis ticks font size
    ax_lw = f_info.ax_lw

    #
    fig, ax = plt.subplots()
    ax.boxplot(dat[:, :, 0], widths = 0.40)

    #
    ax.set_xticklabels(labels)
    ax.set_xlabel('Oscillation Bands', {'fontsize': ax_fs, 'fontweight': 'bold'})
    ax.set_ylabel('Correlation Value', {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Set the top and right side frame & ticks off
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Set linewidth of remaining spines
    ax.spines['left'].set_linewidth(ax_lw)
    ax.spines['bottom'].set_linewidth(ax_lw)

    # Add title to plot
    plt.title('Spatial Analysis', {'fontsize': t_fs, 'fontweight': 'bold'})
