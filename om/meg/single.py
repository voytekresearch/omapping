"""MEG-DATA Analysis Module - Single Subject"""

import os
import csv
import datetime

import numpy as np
import scipy.io as sio
from scipy.stats.stats import pearsonr

from fooof import FOOOFGroup

from om.core.utils import clean_file_list, get_cur_subj, extract_fooof_group
from om.core.errors import DataNotComputedError, UnknownDataSourceError, InconsistentDataError

###################################################################################################
###################################################################################################

class MegSubj(object):
    """Class for a single subject of fooof results for MEG Source PSDs.

    Attributes
    ----------
    dat_source : {'OMEGA', 'HCP'}
        Which database the MEG subject comes from.
    db : OMDB() object
        Database object of paths for the omegamappin project.
    subnum : int
        Subject number of the subject of MEG data that is loaded.
    n_psds : int
        Number of PSDs loaded, corresponds to number of vertices.
    slopes : 1d array
        Slope value for each vertex.
    centers : 2d array
        Oscillation centers, for each vertex [n_verts, n_oscs].
    powers : 2d array
        Oscillation powers, for each vertex [n_verts, n_oscs].
    bws : 2d array
        Oscillation bandwidths, for each vertex [n_verts, n_oscs].
    centers_all : 1d array
        Vector of all oscillation centers, collapsed across vertices.
    powers_all : 1d array
        Vector of all oscillation powers, collapsed across vertices.
    bws_all : 1d array
        Vector of all oscillation bandwidths, collapsed across vertices.
    n_oscs : int
        Number of oscillations extracted for the subject, across all vertices.
    centers_hist : 1d array
        Center frequencies binned into histogram format.
    bands : Osc() object
        Oscillation band definitions.
    oscs : dict
        Oscillations, split up by band, per vertex.
    peaks : dict
        Peak frequencies within each oscillation band.
    comment : str
        Note on the current data, typically subject number, can be used for plotting.
    sex : {'M', 'F'}
        Sex of the current subject.
    age : int
        Age of the current subject.
    osc_count : 1d array
        Number of oscillations extracted per vertex.
    has_data : boolean
        Whether or not data has been loaded to current object.
    has_vertex_oscs : boolean
        Whether data is defined for oscillations at each vertex (not in bands).
    has_vertex_slopes : boolean
        Whether slope values per vertex are loaded.
    has_all_osc : boolean
        Whether data has been flattened into all oscillations.
    has_vertex_bands : boolean
        Whether data has been converted into band specific oscillations, per vertex.
    has_peak_freqs : boolean
        Whether peak frequencies, within oscillatory bands, have been calculated.
    has_demo : boolean
        Whether demographic data has been loaded.
    """

    def __init__(self, db, dat_source, osc=None):
        """Initialize object with omegamappin database.

        Parameters
        ----------
        db : OMDB() object
            Database object for omegamappin project.
        osc : Osc() object, optional
            Oscillation band definitions.
        """

        # Store which db is set
        self.dat_source = dat_source

        # Add database object
        self.db = db

        # Initialize subject number
        self.subnum = int()
        self.n_psds = int()

        # Initialize data arrays
        self.slopes = np.array([])
        self.centers = np.array([])
        self.powers = np.array([])
        self.bws = np.array([])

        # Initialize arrays for all-osc data
        self.centers_all = np.array([])
        self.powers_all = np.array([])
        self.bws_all = np.array([])
        self.n_oscs = np.array([])

        # Initialize list to store all centers histogram
        #  This is used for the subject specific oscillation analysis/plots
        self.centers_hist = []

        # Initialize dictionary to store osc bands
        if osc:
            self.bands = osc.bands
        else:
            self.bands = dict()

        # Initialize dictionary for oscillation bands
        self.oscs = dict()

        # Initialize dictionary for peak frequencies
        self.peaks = dict()

        # Set plot comment
        self.comment = ''

        # Set boolean for what has been run
        self.has_data = False
        self.has_vertex_oscs = False
        self.has_vertex_slopes = False
        self.has_all_osc = False
        self.has_vertex_bands = False
        self.has_peak_freqs = False
        self.has_demo = False

        # Initialize demographic variables
        self.sex = str()
        self.age = int()

        # Initialize oscillation count
        self.osc_count = int()


    def set_bands(self, osc):
        """

        Parameters
        ----------
        osc : ?
            xx
        """

        if self.bands and (self.has_vertex_bands or self.peaks):
            raise InconsistentDataError("Can't change band definitions after they have been used.")
        else:
            self.bands = osc.bands


    def import_fooof(self, subnum, get_demo=True):
        """Import fooof results to MegSubj object.

        Parameters
        ----------
        subnum : int
            Number of the subject to import.
        get_demo : boolean, optional (default = True)
            Whether to load demographic data from csv file.
        """

        # Check if object already has data
        if self.has_data:
            raise InconsistentDataError('Subject object already contains data. Can not add')

        # Set subject number for current data object
        self.subnum = subnum
        self.comment = 'S-' + str(self.subnum)

        # Set up paths, get list of files for available subjects
        files = os.listdir(os.path.join(self.db.fooof_path, self.dat_source))
        files = clean_file_list(files, 'fooof_Vertex')

        # Get specific file name, and set up full file path for specified subject
        cur_subj_file = get_cur_subj(subnum, files)
        cur_subj_path = os.path.join(self.db.fooof_path, self.dat_source, cur_subj_file)

        # Load FOOOF data file
        fg = FOOOFGroup()
        fg.load(cur_subj_path)
        self.centers, self.powers, self.bws, self.slopes, self.n_psds = extract_fooof_group(fg)

        # Update which data is loaded
        self.has_vertex_oscs = True
        self.has_vertex_slopes = True

        # Check how many oscillations per vertex
        self.osc_count = np.zeros([self.n_psds])
        for i in range(self.n_psds):
            self.osc_count[i] = len(np.nonzero(self.centers[i, :])[0])

        # Get demographic data
        if get_demo:
            self.has_demo = True
            self.sex, self.age = _get_demo_csv(self.subnum, self.db.meg_path, self.dat_source)

        # Update boolean to say current subject has data attached
        self.has_data = True


    def all_oscs(self, verbose=True):
        """Flatten osc data to vectors.

        Parameters
        ----------
        verbose : boolean
            Whether to print out status as function runs.

        Notes
        -----
        - When imported, oscillation data is in matrix form [n_vertex, osc_slots].
            This functions converts these matrices into 1-D vectors.
        - This function loses information about which vertex oscillations occur at.
        """

        # Flatten osc data into vectors. Uses C-style row-major order
        self.centers_all = self.centers.flatten('C')
        self.powers_all = self.powers.flatten('C')
        self.bws_all = self.bws.flatten('C')

        # Flattened vectors will have lots of zeros. Get only non-zero indices.
        # TODO: note - currently checks based on powers : hack for dropping 0 amplitudes
        non_zeros = np.nonzero(self.powers_all)
        self.centers_all = self.centers_all[non_zeros]
        self.powers_all = self.powers_all[non_zeros]
        self.bws_all = self.bws_all[non_zeros]

        # Check for nans in BW estimation
        # NOTE: Updated fooof sometimes returns NaN for bw. Check and discard those.
        nans = np.isnan(self.bws_all)
        n_nans = sum(nans)

        # If there are nans, print how many and remove them
        if n_nans > 0:
            if verbose:
                print('Subj:', str(self.subnum), 'Found', str(n_nans), ' NaNs. Removing.')

            # Remove osc's with nan values
            non_nans = np.logical_not(nans)
            self.centers_all = self.centers_all[non_nans]
            self.powers_all = self.powers_all[non_nans]
            self.bws_all = self.bws_all[non_nans]

        # Create the all-centers histogram
        self.centers_hist, _ = np.histogram(self.centers_all, bins=np.arange(3, 40.25, 0.25))

        # Get the number of oscillations
        self.n_oscs = len(self.centers_all)

        # Update boolean that all-oscillations has been computed
        self.has_all_osc = True


    def osc_bands_vertex(self):
        """Groups oscillations at each vertex in distinct frequency bands.

        Stores band specific oscillations in (self.){thetas, alphas, betas, lowgammas}.
        """

        # Check that oscillation bands are defined
        if not self.bands:
            raise DataNotComputedError('Oscillation bands not specified, can not proceed.')

        # Initialize matrices to store oscillations in each band
        for band in self.bands:
            self.oscs[band] = np.zeros([self.n_psds, 4])

        # Loop through each vertex
        for i in range(self.n_psds):

            # Get centers, powers and bws from individual vertex
            centers_temp = self.centers[i, :]
            powers_temp = self.powers[i, :]
            bws_temp = self.bws[i, :]

            # Loop through each oscillation band
            for band in self.bands:

                # Get oscillations in specific band
                self.oscs[band][i, :] = _get_single_osc(centers_temp, powers_temp, bws_temp,
                                                        self.bands[band][0], self.bands[band][1])

        # Update boolean to note that current subject has band specific oscs calculated.
        self.has_vertex_bands = True


    def peak_freq(self, dat, avg='mean'):
        """Calculates the peak frequency for each oscillatory band.

        Parameters
        ----------
        dat : {'all', 'band'}
            Which data to use to calculate peak frequency.
        avg : {'mean', 'median'}, optional
            Which type of averaging to do.
        """

        # Check that oscillation bands are defined
        if not self.bands:
            raise DataNotComputedError('Oscillation bands not specified, can not proceed.')

        # Calculate peak frequency from all oscillations, across all vertices, within band range
        if dat is 'all':

            # Check all osc data has been computed
            if not self.has_all_osc:
                raise DataNotComputedError('All Osc data has not been computed. Can not continue.')

            # Loop through each band, calculating peak frequency
            for band in self.bands:
                self.peaks[band] = _osc_peak_all(
                    self.centers_all, self.bands[band][0], self.bands[band][1], avg)

        # Calculate peak frequency from within band specific maps, which have enforced
        #  a single oscillation per band, per vertex, choosing highest power oscillations
        elif dat is 'band':

            if not self.has_vertex_bands:
                raise DataNotComputedError('Bands not computed per vertex, can not continue.')

            for band in self.bands:

                non_zero_inds = np.nonzero(self.oscs[band][:, 0])

                if avg is 'mean':
                    self.peaks[band] = np.mean(self.oscs[band][non_zero_inds, 0])
                elif avg is 'median':
                    self.peaks[band] = np.median(self.oscs[band][non_zero_inds, 0])

        # Update tracking of what data has been loaded/computed
        self.has_peak_freqs = True


    def calc_osc_param_corrs(self):
        """Calculates correlations between oscillatory parameters.

        Returns
        -------
        corrs_mat : ?
            xx
        ps_mat : ?
            xx
        labels : ?
            xx
        """

        # Check all osc data has been computed
        if not self.has_all_osc:
            raise DataNotComputedError('All Osc data has not been computed. Can not continue.')

        # Set labels for the things being correlated
        labels = ['Centers', 'Powers', 'Bandwidths']

        # Check how many categories there are
        n_label = len(labels)

        # Initialize matrices to store R and p values
        corrs_mat = np.zeros([n_label, n_label])
        ps_mat = np.zeros([n_label, n_label])

        # Calculate correlations between all parameters
        corrs_mat[0, 1], ps_mat[0, 1] = pearsonr(self.centers_all, np.log10(self.powers_all))
        corrs_mat[0, 2], ps_mat[0, 2] = pearsonr(self.centers_all, np.log10(self.bws_all))
        corrs_mat[1, 2], ps_mat[1, 2] = pearsonr(np.log10(self.powers_all),
                                                 np.log10(self.bws_all))

        # Duplicate matrix across the diagonal
        corrs_mat = corrs_mat + corrs_mat.T
        ps_mat = ps_mat + ps_mat.T

        return corrs_mat, ps_mat, labels


    def set_fooof_viz(self):
        """Saves a matfile of freq info to be loaded with Brainstorm for visualization.

        Notes
        -----
        Requires that ... data be
        """

        # Set up paths to save to
        save_name = str(self.subnum) + '_fooof_Viz'
        save_file = os.path.join(self.db.viz_path, save_name)

        # Initialize dictionary, save basic information and slope data
        save_dict = {}
        save_dict['subnum'] = self.subnum
        save_dict['dat_source'] = self.dat_source
        save_dict['save_time'] = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        save_dict['slopes'] = self.slopes

        # Save out oscillation data
        for band in self.bands:
            save_dict[band.lower()] = self.oscs[band]

        # Save the dicionary out to a .mat file
        sio.savemat(save_file, save_dict)


############################################################################################
########################## OMEGAMAPPIN - OM_MD - PUBLIC FUNCTIONS ##########################
############################################################################################

def print_corrs_mat(rs_mat, ps_mat, labels):
    """Prints out correlations, from a given matrix.

    Parameters
    ----------
    rs_mat : 2d array
        Matrix of R values.
    ps_mat : 2d array
        Matrix of p-values.
    labels : list of str
        Labels for the rows & columns.
    """

    # Check how size of the matrix there are
    n_label = len(labels)

    # Loop through the matrix to print out
    for x_dat in range(n_label):
        for y_dat in range(n_label):

            # Skip bottom triangle and diagonal
            if x_dat == y_dat or y_dat < x_dat:
                continue

            # Print out correlation
            print('Corr of ', '{:18}'.format(labels[x_dat]+'-'+labels[y_dat]),
                  ' is ', '{:+1.4f}'.format(rs_mat[x_dat, y_dat]), '    with p-val of ',
                  '{:1.5f}'.format(ps_mat[x_dat, y_dat]))


def print_corrs_vec(rs_vec, ps_vec, labels, desc):
    """Prints out correlations, from a given vector.

    Parameters
    ----------
    rs_vec : 1d array
        Vector of R values.
    ps_vec : 1d array
        Vector of p-values.
    labels : list of str
        Labels for the columns.
    desc : str
        Label for the row.
    """

    # Check the length of the vector
    n_label = len(labels)

    # Loop through vectors, printing out the correlations.
    for cur_label in range(n_label):
        print('Corr of ', '{:20}'.format(labels[cur_label]+'-'+desc), ' is ',
              '{:+1.4f}'.format(rs_vec[cur_label]), '    with p-val of ',
              '{:1.5f}'.format(ps_vec[cur_label]))

#################################################################################################
############################ OMEGAMAPPIN - OM_MD - PRIVATE FUNCTIONS ############################
#################################################################################################

def _get_single_osc(centers, powers, bws, osc_low, osc_high):
    """ Searches for an oscillations of specified frequency band.

    Returns a single oscillation in that band.
    Helper function for osc_per_vertex in MegSubj.

    Parameters
    ----------
    centers : 1d array
        Vector of oscillation centers.
    powers : 1d array
        Vector of oscillation powers.
    bws : 1d array
        Vector of oscillation bandwidths.
    osc_low : float
        Lower bound of frequency band to extract.
    osc_high : float
        Upper bound of frequency band to extract.

    Returns
    -------
    osc_out : array
        Osc data, form - (centers, powers, bws, # oscillations).
    """

    # Find indices of oscillations in the specified range
    osc_inds = (centers >= osc_low) & (centers <= osc_high)

    # Get cen, pow & bw for oscillations in specfied range
    osc_cens = centers[osc_inds]
    osc_pows = powers[osc_inds]
    osc_bws = bws[osc_inds]

    # Get number of oscillations in the frequency band.
    n_oscs = len(osc_cens)

    # Get highest power oscillation in band
    cen, power, bw = _get_single_osc_power(osc_cens, osc_pows, osc_bws)

    return np.array([cen, power, bw, n_oscs])


def _get_single_osc_power(osc_cens, osc_pows, osc_bws):
    """Return the highest power oscillation in a given range.

    Parameters
    ----------
    osc_cens : 1d array
        Vector of oscillation centers.
    osc_pows : 1d array
        Vector of oscillation powers.
    osc_bws : 1d array
        Vector of oscillation bandwidths.

    Returns
    -------
    center : float
        Center frequency value of highest power oscillation.
    power : float
        Power value of highest power oscillation.
    bw : float
        Bandwidth of highest power oscillation.
    """

    # Return zeros if there are no oscillations in given vectors
    if len(osc_cens) == 0:
        return 0., 0., 0.

    # If singular oscillation, return that oscillation
    elif len(osc_cens) == 1:
        return osc_cens, osc_pows, osc_bws

    # If multiple oscillations, return the one with the highest power
    else:
        high_ind = np.argmax(osc_pows)
        return osc_cens[high_ind], osc_pows[high_ind], osc_bws[high_ind]


def _get_demo_csv(subnum, meg_path, dat_source, use_restricted=True):
    """Get demographic information from csv file for specified subject.

    Parameters
    ----------
    subnum : int
        Subject number to get demographic info for.
    meg_path : str
        String for path to csv file.
    dat_source: {'OMEGA', 'HCP'}
        Which database subject is from.

    Returns
    -------
    sex : {'M', 'F'}
        Sex ['M'/'F'] of specified subject.
    age : int
        Age (in whole years) of specified subject.
    """

    # Set up paths for demographic info csv file
    if dat_source is 'OMEGA':
        num_ind = 1
        sex_ind = 4
        age_ind = 7
        res = ''
    elif dat_source is 'HCP':
        if use_restricted:
            num_ind = 0
            sex_ind = 6
            age_ind = 1
            res = '_RESTRICTED_wgender'
        else:
            res = ''
            num_ind = 0
            sex_ind = 3
            age_ind = 4

    else:
        raise UnknownDataSourceError('Unrecognized database source to load from.')

    # Set up name of CSV file to use
    csv_file = os.path.join(meg_path, '00-' + dat_source + '_Subjects' + res + '.csv')

    # Open csv file, loop through looking for right row, grab age & sex information
    with open(csv_file, 'r') as f_name:
        reader = csv.reader(f_name, delimiter=',')

        for row in reader:

            # Find subject line in file
            if row[num_ind] == str(subnum):

                #
                if dat_source is 'OMEGA':
                    sex = row[sex_ind]
                    age = int(row[age_ind])

                elif dat_source is 'HCP' and not use_restricted:
                    sex = row[sex_ind]
                    age_temp = (row[age_ind]).split('-')
                    age = (int(age_temp[0]) + int(age_temp[1]))/2

                else:
                    sex = row[sex_ind]
                    age = int(row[age_ind])

                break

    return sex, age


def _osc_peak_all(centers, osc_low, osc_high, avg='mean'):
    """Find the peak-frequency of a vector of center frequencies.

    Parameters
    ----------
    centers : 1d vector
        Vector of center frequencies to use.
    osc_low : float
        Lower bound of frequency range to check.
    osc_high : float
        Upper bound of frequency range to check.
    avg : {'mean', 'median'}, optional (default = 'mean')
        What kind of average to take.

    Returns
    -------
    peak : float
        Peak frequency value, the average frequency within a given range.
    """

    # Pull out all center frequencies between given range
    osc_inds = (centers > osc_low) & (centers < osc_high)
    osc_cens = centers[osc_inds]

    # Take the average of the center frequencies
    if avg is 'mean':
        peak = np.mean(osc_cens)
    elif avg is 'median':
        peak = np.median(osc_cens)

    return peak
