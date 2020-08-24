"""OM - plots for MEG data."""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from om.plts.fig_info import FigInfo
from om.plts.utils import set_axis_spines, save_figure
from om.core.errors import UnknownDataTypeError

###################################################################################################
###################################################################################################

############################## OM - PLTS - MEGDATA PLOTS ##############################

def plot_exponents(exponents, title, save_out=False, fig_info=FigInfo()):
    """Plots a histogram of the chi values for all vertices.

    Parameters
    ----------
    exponents : 1d array
        A vector of aperiodic exponent values to plot.
    title : str
        A string to append to the title.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Plot Settings
    n_bins = 150             # Number of bins for histograms
    t_fs = fig_info.t_fs           # Title font size
    ax_fs = fig_info.ax_fs         # Axis label font size
    ti_fs = fig_info.ti_fs         # Axis ticks font size

    # Set up plot
    fig, ax = plt.subplots(figsize=[6, 4])

    # Create histogram
    ax.hist(exponents, n_bins, color='#40425e')

    # Add title
    if fig_info.add_title:
        plt.title('Exponents - ' + title, {'fontsize': t_fs})

    # Add axis labels
    ax.set_xlabel('Exponent', {'fontsize': ax_fs})
    ax.set_ylabel('Count', {'fontsize': ax_fs})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Set spines
    set_axis_spines(ax, lw=fig_info.ax_lw)

    # Set x-lims
    ax.set_xlim(0.0, 2.0)

    save_figure(save_out, '101-' + title + '_Exponents', fig_info)


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
    fig_info = FigInfo()

    # Plot Settings
    n_bins = 25              # Number of bins for histograms
    t_fs = fig_info.t_fs           # Title font size
    ax_fs = fig_info.ax_fs         # Axis label font size
    ti_fs = fig_info.ti_fs         # Axis ticks font size

    # Create histogram
    plt.hist(osc_count, n_bins, range=[0, 8])

    # Add title
    if fig_info.add_title:
        plt.title('# Oscillations per Vertex', {'fontsize': t_fs, 'fontweight': 'bold'})

    # Add axis labels
    plt.xlabel('# Oscillations', {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    save_figure(save_out, '102-OscCount', fig_info)


def plot_all_oscs(centers_all, powers_all, bws_all, title,
                  save_out=False, fig_info=FigInfo()):
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

    # Plot Settings
    n_bins = 160                 # Number of bins for histograms
    st_fs = fig_info.t_fs          # Super Title Font Size
    sp_fs = fig_info.sp_fs         # Subplot Title Font Size
    ax_fs = fig_info.ax_fs         # Axis Label Font Size
    ti_fs = fig_info.ti_fs         # Axis ticks font size

    # Set up subplots
    fig, ax = plt.subplots(3, 1, figsize=(15, 15))

    # Set plot super-title
    if fig_info.add_title:
        plt.suptitle('Distributions of Oscillatory Parameters - ' + title,
                     fontsize=st_fs, fontweight='bold')

    # Subplot 1 - Center Frequency
    ax[0].hist(centers_all, n_bins)
    ax[0].set_title('Center Frequency', {'fontsize': sp_fs})
    ax[0].set_xlabel('Frequency', {'fontsize': ax_fs})
    ax[0].set_ylabel('Count', {'fontsize': ax_fs})

    # Subplot 2 - Power
    ax[1].hist(np.log10(powers_all), n_bins)
    ax[1].set_title('Oscillatory Power', {'fontsize': sp_fs})
    ax[1].set_xlabel('Log Power', {'fontsize': ax_fs})
    ax[1].set_ylabel('Count', {'fontsize': ax_fs})

    # Subplot 3 - Bandwidth
    ax[2].hist(bws_all, n_bins)
    ax[2].set_title('Band Width', {'fontsize': sp_fs})
    ax[2].set_xlabel('Bandwidth (Hz)', {'fontsize': ax_fs})
    ax[2].set_ylabel('Count', {'fontsize': ax_fs})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Adjust subplot spacing
    plt.subplots_adjust(hspace=0.4)

    save_figure(save_out, '103-AllOscs', fig_info)


def plot_all_oscs_single(data, data_type, title, n_bins=160,
                         figsize=(15, 5), save_out=False, fig_info=FigInfo()):
    """Create a plot for a single oscillation parameter.

    Parameters
    ----------
    data : 1d array
        Vector of oscillation parameter data.
    data_type : {0, 1, 2}
        Int refers to which osc parameter is being plotted.
            Key: {0:'Center Frequency', 1:'Power', 2:'Bandwidth'
    title : str
        A string to append to the title.
    n_bins : int, optional (default = 160)
        Number of bins to use for the plot.
    figsize : tuple, optional (default = (15, 5))
        Size of the figure to make.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Plot Settings
    t_fs = fig_info.t_fs           # Super Title Font Size
    ax_fs = fig_info.ax_fs         # Axis Label Font Size
    ti_fs = fig_info.ti_fs         # Axis ticks font size

    # Set up for which data type
    if data_type is 0:
        data_title = 'Center Frequency'
        xlab = 'Frequency'
    elif data_type is 1:
        data_title = 'Power'
        xlab = 'Log Power'
        data = np.log10(data)
    elif data_type is 2:
        data_title = 'Bandwidth'
        xlab = 'Bandwidth (Hz)'
    else:
        raise UnknownDataTypeError('Data type not understood.')

    # Set up plot
    fig, ax = plt.subplots(figsize=figsize)

    # Subplot 1 - Center Frequency
    ax.hist(data, n_bins, color='#40425e', alpha=0.9)

    # Set the title
    if fig_info.add_title:
        plt.title(data_title + ' ' + title, {'fontsize': t_fs, 'fontweight': 'bold'})

    # Add axis labels
    ax.set_xlabel(xlab, {'fontsize': ax_fs})
    ax.set_ylabel('Count', {'fontsize': ax_fs})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Set spines
    set_axis_spines(ax, lw=fig_info.ax_lw)

    # Hard code x-lims
    #ax.set_xlim(3, max(data)+0.2)
    ax.set_xlim(2.5, 32)

    save_figure(save_out, '104-' + data_title + '_AllOscs', fig_info)


def plot_osc_param_comparison(centers_all, powers_all, bws_all, title,
                              save_out=False, fig_info=FigInfo()):
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

    # Plot Settings
    st_fs = fig_info.t_fs              # Super Title Font Size
    sp_fs = fig_info.sp_fs             # Subplot Title Font Size
    ax_fs = fig_info.ax_fs             # Axis Label Font Size
    ti_fs = fig_info.ti_fs             # Axis ticks font size
    vis_opac = 0.1                   # Alpha value for plotted data

    # Set up subplots
    fig, ax = plt.subplots(3, 1, figsize=(15, 15))

    # Set plot super-title
    if fig_info.add_title:
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

    save_figure(save_out, '105-OscComparison', fig_info)


def plot_corr_matrix(corr_dat, labels, save_out=False, fig_info=FigInfo()):
    """Plot correlation data.

    Parameters
    ----------
    corr_data : 2d array
        Matrix of correlation data to plot.
    labels : list of str
    	Labels for the rows & columns of `corr_data`.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Plot Settings
    t_fs = fig_info.t_fs              # Title Font Size
    ti_fs = fig_info.ti_fs            # Axis ticks font size

    # Set colormap to use
    cmap = plt.get_cmap('seismic')

    # Create the plot
    #im = plt.matshow(corr_data, vmin=-1, vmax=1, cmap=cmap, interpolation='none')
    im = plt.matshow(corr_data, vmin=-1, vmax=1, cmap=cmap, interpolation='nearest')
    # Notes on using nearest here:
    #   https://github.com/matplotlib/matplotlib/issues/2972/

    # Add title
    if fig_info.add_title:
        plt.title('Osc Band Correlations', {'fontsize': t_fs, 'fontweight': 'bold'}, y=1.15)

    # Set tick labels
    nums = list(range(len(labels)))
    plt.xticks(nums, labels, rotation=45, ha='left')
    plt.yticks(nums, labels)

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Add a colorbar - add padding to offset further from plot
    plt.colorbar(pad=0.15)

    save_figure(save_out, '106-CorrMat', fig_info)


def plot_corr_matrix_tri(corr_data, labels, save_out=False):
    """   """

    # Generate a mask for the upper triangle
    mask = np.zeros_like(corr_data, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    # Generate a custom diverging colormap
    cmap = sns.color_palette("coolwarm", 7)

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corr_data, mask=mask, cmap=cmap, annot=True, square=True,
                vmin=-1, vmax=1, xticklabels=labels, yticklabels=labels)

    save_figure(save_out, 'XX-corr_mat_tri', fig_info)


def plot_peak_boxplot(peaks, osc_band, save_out=False, fig_info=FigInfo()):
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

    # Plot Settings
    t_fs = fig_info.sp_fs              # Title Font Size
    ti_fs = fig_info.ti_fs            # Axis ticks font size

    # Initialize the figure
    fig, ax = plt.subplots(figsize=(2, 5))

    # Make the boxplot
    ax.boxplot(peaks, widths=0.45)

    # Set spines
    set_axis_spines(ax, lw=fig_info.ax_lw)

    # Add a title
    #if fi.add_title:
    plt.title(osc_band, fontsize=t_fs, fontweight='bold', y=1.08)

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Remove x ticks
    plt.xticks([])

    save_figure(save_out, '107-' + osc_band + '_boxplot', fig_info)


def plot_peak_boxplot_all(peaks, bands, save_out=False, fig_info=FigInfo()):
    """Plot boxplots of peak frequencies for all oscillation bands.

    Parameters
    ----------
    peaks : ?
        xx
    bands : ?
        xx

    TODO: Figure out how to space out subplots a little
    """

    # Plot Settings
    t_fs = fig_info.sp_fs             # Title Font Size
    ti_fs = fig_info.ti_fs            # Axis ticks font size
    ax_lw = fig_info.ax_lw            # Axis line weight

    n_bands = len(bands)
    fig, axes = plt.subplots(1, n_bands, figsize=(10, 6))

    for ind, band in enumerate(bands):
        ax[ind].boxplot(peaks[band], widths=0.40)

        # Add title and set y lims using bands and data
        axes[ind].set_title(band, fontsize=t_fs, fontweight='bold', y=1.04)
        axes[ind].set_ylim([peaks[band].min()-0.2, peaks[band].max()+0.2])

    for ax in axes:

        # Set spines & tick params
        set_axis_spines(ax, lw=fig_info.ax_lw)
        ax.tick_params(axis='both', which='major', labelsize=ti_fs)
        ax.set_xticks([])

    save_figure(save_out, 'XXX-PeakFreqs_All', fig_info)


def plot_freq_corr(fs, corr_vec, p_vec, save_out=False, fig_info=FigInfo()):
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

    # Plot settings
    t_fs = fig_info.t_fs
    ax_fs = fig_info.ax_fs
    ti_fs = fig_info.ti_fs         # Axis ticks font size
    ax_lw = fig_info.ax_lw

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
    if fig_info.add_title:
        plt.suptitle('Correlation Adjacent Frequency Bands', fontsize=t_fs, fontweight='bold')

    # Set axis limits
    ax.set_xlim([0, 30.5])
    ax.set_ylim([-1.0, 1.0])

    # Set axis labels
    ax.set_xlabel('Frequency', fontsize=ax_fs, fontweight='bold')
    ax.set_ylabel('R Value', fontsize=ax_fs, fontweight='bold')

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Set spines
    set_axis_spines(ax, lw=fig_info.ax_lw)

    save_figure(save_out, '108-FreqCorr', fig_info)


def plot_age_peak(age, peak_theta, peak_alpha, peak_beta,
                  save_out=False, fig_info=FigInfo()):
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

    # Plot settings
    st_fs = fig_info.t_fs
    sp_fs = fig_info.sp_fs
    ti_fs = fig_info.ti_fs         # Axis ticks font size

    # Set up subplots
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))

    # Set plot super-title
    if fig_info.add_title:
        plt.suptitle('Peak Frequency / Age Comparisons', fontsize=st_fs, fontweight='bold')

    # Theta
    ax[0].plot(age, peak_theta, '.')
    ax[0].set_title('Theta', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[0].tick_params(axis='both', which='major', labelsize=ti_fs)

    # Alpha
    ax[1].plot(age, peak_alpha, '.')
    ax[1].set_title('Alpha', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[1].tick_params(axis='both', which='major', labelsize=ti_fs)

    # Beta
    ax[2].plot(age, peak_beta, '.')
    ax[2].set_title('Beta', {'fontsize': sp_fs, 'fontweight': 'bold'})
    ax[2].tick_params(axis='both', which='major', labelsize=ti_fs)

    # Set ticks font size
    #plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    save_figure(save_out, '109-AgePeak', fig_info)


def plot_scatter(d1, d2, labels=[None, None], title=None,
                 save_out=False, fig_info=FigInfo()):
    """Create a scatter plot.

    Parameters
    ----------
    d1 : 1d array
        Vector of x-axis data.
    d2 : 1d array
        Vector of y-axis data.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Plot settings
    t_fs = fig_info.t_fs
    ti_fs = fig_info.ti_fs
    ax_fs = fig_info.ax_fs

    # Make the plot
    plt.plot(d1, d2, '.', markersize=10, alpha=0.8)

    # Add title
    if fig_info.add_title:
        plt.title(title, fontsize=t_fs, fontweight='bold')

    # Add axis labels
    plt.xlabel(labels[0], {'fontsize': ax_fs, 'fontweight': 'bold'})
    plt.ylabel(labels[1], {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    save_figure(save_out, '222' + title, fig_info)


def plot_osc_profiles(centers_hist, n_subj='all',
                       figsize=(15, 5), save_out=False, fig_info=FigInfo()):
    """Creates a plot showing all subjects oscillation profiles.

    Parameters
    ----------
    centers_hist : list of list
        Contains the oscillation profile for each subject.
    n_subj : 'all' or int
        If all: plots all subj. If int, plots n_sing random individual subjects.
    save_out : boolean, optional (default = False)
        Whether to save out a copy of the figure.
    """

    # Plot settings
    t_fs = fig_info.t_fs
    ax_fs = fig_info.ax_fs
    ti_fs = fig_info.ti_fs         # Axis ticks font size
    ax_lw = fig_info.ax_lw

    # xx
    ind_lw = 1.0
    avg_lw = 3.5
    alpha = 0.80

    # Initialize plot
    fig, ax = plt.subplots(figsize=figsize)

    # Initialize a frequency vector for the x-axis
    freqs = np.arange(3.125, 40.125, 0.25)

    # Set inds to plot
    if n_subj is 'all':
        subj_inds = np.range(len(centers_hist))
    else:
        subj_inds = np.random.choice(range(len(centers_hist)), n_subj, replace=False)

    # Loop through all subjects, adding profile to plot
    for ind, hist in enumerate(centers_hist):
        if ind in subj_inds:
            ax.plot(freqs, hist, linewidth=ind_lw, alpha=alpha)

    # Add the average profile
    ax.plot(freqs, np.median(centers_hist, 0), 'k', linewidth=avg_lw)

    # Add title
    if fig_info.add_title:
        plt.title('Individual Oscillation Profiles', {'fontsize': t_fs, 'fontweight': 'bold'})

    # Set spines
    set_axis_spines(ax, lw=fig_info.ax_lw)

    # Add axis labels
    ax.set_xlabel('Frequency', {'fontsize': ax_fs})
    ax.set_ylabel('Count', {'fontsize': ax_fs})

    # Set ticks font size
    plt.tick_params(axis='both', which='major', labelsize=ti_fs)

    # Set the x axis limits
    ax.set_xlim([3, 30])

    save_figure(save_out, '111-OscillationProfiles', fig_info)


def plot_space_comp(oscs, verts, band, subj=0, osc_param=0, space_param=1,
                    save_out=False, fig_info=FigInfo()):
    """

    TODO: Fill in.

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

    # Plot settings
    t_fs = fig_info.t_fs
    ax_fs = fig_info.ax_fs

    #
    space = verts[:, space_param]
    sort_inds = np.argsort(space)
    freqs = [data if data > 0 else None for data in oscs[band][sort_inds, osc_param, subj]]

    fig = plt.figure()

    plt.plot(space, freqs, '.', ms=3.5, alpha=0.75)

    plt.xlim([space.min()-0.05, space.max()+0.05])

    plt.title(band, {'fontsize': t_fs})

    # NOTE: These are wrong / need to change with osc & space params
    #plt.xlabel('Posterior -> Anterior', {'fontsize': ax_fs, 'fontweight': 'bold'})
    #plt.ylabel('Center Frequency', {'fontsize': ax_fs, 'fontweight': 'bold'})

    save_figure(save_out, '1XX-OscillationSpaceSingleSubject', fig_info)


def plot_space_comp_all():
    """   """
    pass


def plot_osc_space_corr_boxplot(data, labels, save_out=False, fig_info=FigInfo()):
    """

    Parameters
    ----------
    data : ?
        xx
    labels : ?
        xx
    """

    # Plot settings
    t_fs = fig_info.t_fs
    ax_fs = fig_info.ax_fs
    ti_fs = fig_info.ti_fs         # Axis ticks font size

    fig, ax = plt.subplots()
    ax.boxplot(data[:, :, 0], widths = 0.40)

    ax.set_xticklabels(labels, {'fontsize': 16})
    #ax.set_xlabel('Oscillation Bands', {'fontsize': ax_fs, 'fontweight': 'bold'})
    #ax.set_ylabel('Correlation Value', {'fontsize': ax_fs, 'fontweight': 'bold'})

    # Set spines
    set_axis_spines(ax, lw=fig_info.ax_lw)

    #plt.title('Spatial Analysis', {'fontsize': t_fs, 'fontweight': 'bold'})

    save_figure(save_out, '1XX-OscillationSpaceCorrs', fig_info)
