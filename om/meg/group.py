"""MEG-DATA Analysis Module - Group"""

import os
import pickle
import datetime

import numpy as np
import scipy.io as sio
from scipy.stats.stats import pearsonr

from om.meg.single import MegSubj
from om.core.osc import check_bands
from om.core.errors import DataNotComputedError, InconsistentDataError, UnknownDataTypeError

###################################################################################################
###################################################################################################

class MegGroup(MegSubj):
    """A class to store OMEGA data from multiple subjects.

    Holds all oscillations, regardless of spatial location.

    Attributes
    ----------
    n_subjs : int
        The number of subjects included in the group data.
    subjs : list of int
        List of the subject numbers included in current group data.
    bands : Osc() object
        Stores labels and band definitions of oscillation bands.
    n_oscs_tot : int
        Total number of oscillations found across all subjects.
    comment : str
        A note about the data, label about data that is loaded.
    gr_oscs : dict
        All oscillations, in bands, for all subjects [n_verts, n_oscs, n_subjs].
    osc_probs : dict
        Oscillation probability for each oscillation band, for each vertex.
    osc_pow_ratios : dict
        Oscillation power ratios for each oscillation band, for each vertex.
    osc_scores : dict
        Oscillation scores for each oscillation band, for each vertex.
    vert_exponents : 2d array
        Aperiodic exponent values for each subject, at each vertex [n_verts, n_subjs].
    exponent_gr_avg : 1d array
        Average aperiodic exponent values across subjects for each vertex.
    osc_prob_done : boolean
        Whether oscillation probability has been calculated.
    osc_power_done : boolean
        Whether oscillation power ratio has been calculated.
    osc_score_done : boolean
        Whether oscillation score has been calculated.
    """

    def __init__(self, db, osc):
        """Initialize object with omegamappin database, and oscillation definitions.

        Parameters
        ----------
        db : OMDB() object
            Database object for omegamappin project.
        osc : Osc() object
            Object to store oscillatory band definitions.
        """

        # Initialize from MegSubj() object
        MegSubj.__init__(self, db, 'both')

        # Initialize groups subject variables
        self.n_subjs = int()
        self.subjs = []

        # Update variable types for demographic
        self.sex = list()
        self.age = np.array([])

        # Set definition of oscillation bands used for the group
        self.bands = osc.bands

        # Initialize count of total oscillations, across all subjects
        self.n_oscs_tot = int()

        # Set comment for data, can be used for plotting
        self.comment = 'Group'

        # Initialize dictionary for oscillation band data
        self.gr_oscs = dict()

        # Initilaize dictionary to store oscillation probabilities
        self.osc_probs = dict()

        # Initialize dict to store oscillation power ratios
        self.osc_pow_ratios = dict()

        # Initialize to store oscillation scores
        self.osc_scores = dict()

        # Initialize vars to store exponent values
        self.vert_exponents = np.array([])
        self.exponent_gr_avg = np.array([])

        # Set booleans for what has been run
        self.osc_prob_done = False
        self.osc_power_done = False
        self.osc_score_done = False


    def __len__(self):

        return self.n_subjs


    def add_subject(self, new_subj, add_vertex_oscs=False, add_vertex_exponents=False,
                    add_all_oscs=False, add_vertex_bands=False, add_peak_freqs=False,
                    add_demo=False):
        """Adds a new subject to the MegGroup object.

        Parameters
        ----------
        new_subj : MegSubj() Object
            MEG subject (instance of MegSubj)
        add_vertex_oscs : boolean, optional (default: False)
            Whether to add all oscillations, across vertices.
        add_vertex_exponents : boolean, optional (default: False)
            Whether to add the aperiodic exponents.
        add_all_oscs : boolean, optional (default: False)
            Whether to add the vectors of all oscillations, collapsed across vertices.
        add_vertex_bands : boolean, optional (default: False)
            Whether to add the oscillation band data, across vertices.
        add_peak_freqs : boolean, optional (default: False)
            Whether to add peak frequencies.
        add_demo : boolean, optional (default: False)
            Whether to add demographic information.
        """

        # Check if subject has data
        if not new_subj.has_data:
            raise DataNotComputedError("Empty meg data object. Cannot add data.")

        # Add oscillations per vertex
        if add_vertex_oscs:

            # Check new subject has relevant data
            if not new_subj.has_vertex_oscs:
                raise DataNotComputedError('New subject does not have vertex osc data.')

            if not self.has_data:

                # Add data to group object
                self.centers = new_subj.centers
                self.powers = new_subj.powers
                self.bws = new_subj.bws

                # Update that group contains this data
                self.has_vertex_oscs = True

            else:

                # Check that group has data defined
                if not self.has_vertex_oscs:
                    raise DataNotComputedError('MEG Group does not include vertex osc data.')

                # Add data to group object
                self.centers = np.dstack([self.centers, new_subj.centers])
                self.powers = np.dstack([self.powers, new_subj.powers])
                self.bws = np.dstack([self.bws, new_subj.bws])

        # Add exponents per vertex
        if add_vertex_exponents:

            # Check new subject has relevant data
            if not new_subj.has_vertex_exponents:
                raise DataNotComputedError('New subject does not have vertex exponent data.')

            if not self.has_data:

                # Add data to group object
                self.vert_exponents = new_subj.exponents

                # Update that group contains this data
                self.has_vertex_exponents = True

            else:
                # Check that group has data defined
                if not self.has_vertex_exponents:
                    raise DataNotComputedError('MEG Group does not include vertex exponent data.')

                # Add data to group object
                self.vert_exponents = np.hstack([self.vert_exponents, new_subj.exponents])

        # Add All-Osc Data
        if add_all_oscs:

            # Check that new subject has all_osc data available
            if not new_subj.has_all_osc:
                raise DataNotComputedError('New subject does not have all osc data.')

            # Check that group has data defined
            if self.has_data:
                if not self.has_all_osc:
                    raise DataNotComputedError('MEG Group does not include all osc data.')

            # Add oscillation parameters to current data
            self.centers_all = np.append(self.centers_all, new_subj.centers_all)
            self.bws_all = np.append(self.bws_all, new_subj.bws_all)
            self.powers_all = np.append(self.powers_all, new_subj.powers_all)
            self.exponents = np.append(self.exponents, new_subj.exponents)

            # Add centers hist
            self.centers_hist.append(new_subj.centers_hist)

            # Update count of total number of oscillations
            self.n_oscs = np.append(self.n_oscs, new_subj.n_oscs)
            self.n_oscs_tot = len(self.centers_all)

            # If first subject, update what kind of data is loaded
            if not self.has_data:
                self.has_all_osc = True

        # Add band-specific data
        if add_vertex_bands:

            # Check that new subject has vertex bands data
            if not new_subj.has_vertex_bands:
                raise DataNotComputedError('New subject does not have vertex band data.')

            # Check that new subject has same bands defined
            _ = check_bands([self.bands, new_subj.bands])

            # Add new subject to group oscillations
            if not self.has_data:

                # Add data to group object
                self.gr_oscs = new_subj.oscs

                # Update that group contains this data
                self.has_vertex_bands = True

            else:

                # Check that group has data defined
                if not self.has_vertex_bands:
                    raise DataNotComputedError('MEG Group does not include vertex band data.')

                # Add data to group object
                for band in self.bands:
                    self.gr_oscs[band] = np.dstack([self.gr_oscs[band], new_subj.oscs[band]])

        # Add oscillation peak data
        if add_peak_freqs:

            # Check that new subject has peak frequency data
            if not new_subj.has_peak_freqs:
                raise DataNotComputedError('New subject does not have peak freq data.')

            # Check that new subject has same bands defined
            _ = check_bands([self.bands, new_subj.bands])

            # Add new subject to peak frequencies
            if not self.has_data:

                # Add data to group object
                self.peaks = new_subj.peaks

                # Update that group contains this data
                self.has_peak_freqs = True

            else:

                # Check that group has data defined
                if not self.has_peak_freqs:
                    raise DataNotComputedError('MEG Group does not include peak freq data.')

                # Add data to group object
                for band in self.bands:
                    self.peaks[band] = np.append(self.peaks[band], new_subj.peaks[band])

        # Add demographic data
        if add_demo:

            # Check that incoming subject has demo data
            if not new_subj.has_demo:
                raise DataNotComputedError('Demographic data not available')

            # Check that group has data defined
            if self.has_data:
                if not self.has_demo:
                    raise DataNotComputedError('MEG Group does not include demo data.')

            # Add demographic data to group object
            self.sex.append(new_subj.sex)
            self.age = np.append(self.age, new_subj.age)

            # If first subject, update what kind of data is loaded
            if not self.has_data:
                self.has_demo = True

        # If first subject, update that object has data
        if self.n_subjs == 0:
            self.has_data = True

        # Update subj count and subject number list
        self.n_subjs += 1
        self.subjs.append(new_subj.subnum)

        # Check consistency of group data
        self.check_consistency()


    def check_consistency(self):
        """Check for consistency of data loaded in group object."""

        n_vertices = 7501

        if self.n_subjs != len(self.subjs):
            raise InconsistentDataError('Discrepancy in subject numbers.')

        if self.has_vertex_oscs:
            if self.n_subjs > 1:
                assert self.centers.shape == (n_vertices, 8, self.n_subjs)
                assert self.powers.shape == (n_vertices, 8,self.n_subjs)
                assert self.bws.shape == (n_vertices, 8,self.n_subjs)

        if self.has_vertex_exponents:
            assert self.vert_exponents.shape == (n_vertices, self.n_subjs)

        if self.has_all_osc:
            pass

        if self.has_vertex_bands:
            pass

        if self.has_peak_freqs:
            pass

        if self.has_demo:
            pass


    def group_exponent(self, avg='mean'):
        """Calculates the average exponent value for each vertex, across subjects.

        Parameters
        ----------
        avg : {'mean', 'median'}, optional
            How to average across the group.
        """

        # Calculate the average exponent value per vertex
        if avg is 'mean':
            self.exponent_gr_avg = np.mean(self.vert_exponents, 1)
        elif avg is 'median':
            self.exponent_gr_avg = np.median(self.vert_exponents, 1)


    def osc_prob(self):
        """Calculates the probability of an osc in a specific band.

         This is done per vertex, across subjects.
         """

        # Check if vertex data is set
        if not self.has_vertex_bands:
            raise DataNotComputedError('Vertex oscillation bands data not available.')

        # For each oscillation band, compute the probability of an oscillation in that band - NEW
        for band in self.bands:
            self.osc_probs[band] = _osc_prob(self.gr_oscs[band])

        # Update boolean that oscillation probability has been computed
        self.osc_prob_done = True


    def osc_power(self):
        """Calculate the oscillation power ratio for each frequency band."""

        # Check if vertex data is set
        if not self.has_vertex_bands:
            raise DataNotComputedError('Vertex oscillation bands data not available.')

        # Compute power ratio for each oscillation band
        for band in self.bands:
            self.osc_pow_ratios[band] = _osc_pow_ratio(self.gr_oscs[band])

        # Set boolean that oscillation score has been computed.
        self.osc_power_done = True


    def osc_score(self):
        """Calculate the oscillation score for each frequency band.

        The oscillation score is ....
        """

        # Check if oscillation probability & power ratios are calculated.
        #  Can not proceed if they are not.
        _ = self._get_map_type('prob')
        _ = self._get_map_type('power')

        # Compute oscillation score for each oscillation band
        for band in self.bands:
            self.osc_scores[band] = self.osc_pow_ratios[band] * self.osc_probs[band]

        # Set boolean that oscillation score has been computed.
        self.osc_score_done = True


    def osc_map_corrs(self, map_type):
        """Calculates the correlations between oscillation probabilities or scores.

        Parameters
        ----------
        map_type : {'prob', 'score', 'power'}
            Which map data type to save out.

        Returns
        -------
        corrs_mat : 2d array
            Correlation R-values matrix, across all oscillation bands.
        ps_mat : 2d array
            Correlations p-values matrix, across all oscillation bands.
        sorted_bands : list of str
            Oscillation band labels, sorted into order.
        """

        # Check how many oscillation bands are defined
        n_bands = len(self.bands)

        # Initialize matrices to store correlation results
        corrs_mat = np.zeros([n_bands, n_bands])
        ps_mat = np.zeros([n_bands, n_bands])

        # Get oscillation bands in order
        sorted_bands, sort_inds = _band_sort(self.bands)

        # Set which map to run
        dat = self._get_map_type(map_type)

        # Loop through all bands, computing correlations between them
        for i in range(n_bands):
            for j in range(n_bands):
                corrs_mat[i, j], ps_mat[i, j] = pearsonr(
                    dat[sorted_bands[i]],
                    dat[sorted_bands[j]])

        # Set diagonals to zero - where band is correlated with itself
        np.fill_diagonal(corrs_mat, 0)
        np.fill_diagonal(ps_mat, 0)

        return corrs_mat, ps_mat, sorted_bands


    def calc_osc_peak_age(self):
        """Compares age and peak frequency within frequency bands.

        NOTE: ADD CHECKS THAT REQUIRED DATA HAS BEEN COMPUTED.

        Returns
        -------
        corrs_mat : 1d array
            Correlations R-values comparing age to oscillations.
        ps_mat : 1d array
            Correlations p-values from comparing age to oscillations.
        sorted_bands : list of str
            Oscillation band labels, sorted into order.
        """

        # Check how many bands there are
        n_bands = len(self.bands)

        # Initialize matrices to store correlation results
        corrs_mat = np.zeros([n_bands])
        ps_mat = np.zeros([n_bands])

        # Get oscillation bands in order
        sorted_bands, sort_inds = _band_sort(self.bands)

        # Loop through all oscillation peaks, calculation correlation with age
        for i in range(n_bands):
            corrs_mat[i], ps_mat[i] = pearsonr(
                self.age, self.peaks[sorted_bands[sort_inds[i]]])

        return corrs_mat, ps_mat, sorted_bands


    def save_gr_exponent(self, file_name):
        """Saves out the average group exponent results.

        Parameters
        ----------
        file_name : str
            File name to save group exponent file as.
        """

        # Set up
        pickle_file_name = file_name + '.p'
        pickle_save_name = os.path.join(self.db.maps_path, 'Exponents', pickle_file_name)

        # Check current time for when file is saved
        cur_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # Collect data together to save out
        dat_out = dict({'dat_source': self.dat_source,
                        'exponents': self.exponent_gr_avg,
                        'n_subjs': self.n_subjs,
                        'save_time': cur_time})

        # Save out with pickle
        pickle.dump(dat_out, open(pickle_save_name, 'wb'))


    def save_map(self, map_type, file_name):
        """Save oscillation map data out to disc.

        Parameters
        ----------
        map_type : {'prob', 'score', 'power'}
            Which map data type to save out.
        file_name : str
            String to add to the file name.
        """

        # Set which map to run
        dat = self._get_map_type(map_type)

        # Check current time for when file is saved
        cur_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # Set file name, and create full file path
        pickle_file_name = file_name + '_Osc_' + map_type + '.p'
        pickle_save_name = os.path.join(self.db.maps_path, 'Oscs', pickle_file_name)

        # Collect data together to save out
        dat_out = dict({'dat_source': self.dat_source,
                        'map_type': map_type,
                        'osc_dat': dat,
                        'bands':self.bands,
                        'n_subjs': self.n_subjs,
                        'save_time': cur_time})

        # Save out with pickle
        pickle.dump(dat_out, open(pickle_save_name, 'wb'))


    def set_exponent_viz(self):
        """Saves out a matfile, of the group average exponent, for visualization."""

        # Set up paths to save to
        save_name = 'Group_Exponents'
        save_file = os.path.join(self.db.viz_path, save_name)

        # Save desired outputs into a dictionary
        save_dict = {}
        save_dict['exponents'] = self.exponent_gr_avg
        save_dict['dat_source'] = self.dat_source
        save_dict['n_subjs'] = self.n_subjs
        save_dict['save_time'] = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # Save the dicionary out to a .mat file
        sio.savemat(save_file, save_dict)


    def set_fooof_viz(self):
        """Set FOOOF features to visualize.

        TODO
        """

        pass


    def set_map_viz(self, map_type, file_name):
        """Save out an oscillation map for visualization with Brainstorm.

        Parameters
        ----------
        map_type : {'prob', 'score', 'power'}
            Which map data type to set as viz.
        file_name : str
            Label to attach to file name to be saved out.
        """

        # Set which map to run
        dat = self._get_map_type(map_type)

        # Set up the save name
        save_name = file_name + '_group_osc_' + map_type + '_viz'

        # Set up paths to save to
        save_file = os.path.join(self.db.viz_path, save_name)

        # Initialize dictionary to save out, and save basic info
        save_dict = {}
        save_dict['dat_source'] = self.dat_source
        save_dict['map_type'] = map_type
        save_dict['n_subjs'] = self.n_subjs
        save_dict['save_time'] = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # Add maps to save dictionary
        for band in self.bands:
            save_dict[band.lower() + '_' + map_type] = dat[band]

        # Save out the dictionary to a mat file
        sio.savemat(save_file, save_dict)

    def _get_map_type(self, map_type):
        """Pull out specific MEG map type.

        Parameters
        ----------
        map_type : {'prob', 'score', 'power'}
            Oscillation map type to pull out.
        """

        # Check if requested map is prob, and if it is calculated
        if map_type is 'prob':

            # Check if oscillation probabilities have been calculated.
            if not self.osc_prob_done:
                raise DataNotComputedError("Oscillation probability not computed - can not proceed.")

            dat = self.osc_probs

        # Check if requested map is score, and if it is calculated
        elif map_type is 'score':

            # Check if oscillation score has been calculated.
            if not self.osc_score_done:
                raise DataNotComputedError("Oscillation score not computed - can not proceed.")

            dat = self.osc_scores

        # Check if requested map is power ratio, and if it is calculated
        elif map_type is 'power':

            # Check if oscillation power map has been calculated.
            if not self.osc_power_done:
                raise DataNotComputedError("Oscillation power map not computed - can not proceed.")

            dat = self.osc_pow_ratios

        # Raise an error if requested type doensn't match a known map type
        else:
            raise UnknownDataTypeError('Map type not understood.')

        return dat

#########################################################################################################
################################## OMEGAMAPPIN - MEG GROUP - FUNCTIONS ##################################
#########################################################################################################

def freq_corr_group(centers, f_win, f_step=1):
    """Calculates the correlation between adjacent frequency bands.

    Parameters
    ----------
    centers : 3d array [n_verts, n_slots, n_subjs]
        Center frequencies of oscillations across all vertices & subjects.
    f_win : float
        Size of frequency window to use.
    f_step : float
        Increment to step by.

    Returns
    -------
    corr_vec : 1d array
        Vector of the correlation coefficients between all adjacent frequency bands.
    p_vec : 1d array
        Vector of the p-values for the correlations between adjacent frequency bands.
    freqs : 1d array
        Vector of frequencies of the correlations (each value is first frequency of first bin).
            Each frequency 'f' reflects the correlation of [f:f+f_win, f+f_win:f+2*f_win].
    """

    # Get # vertices, # of subjects to loop through
    ind_step = int(f_win / f_step)
    [n_vertex, n_slots, n_subj] = np.shape(centers)

    # Initialize variables for freqs, # of freqs, and matrix to store probability
    freqs = np.arange(3, 40-f_win, f_step)
    n_freqs = len(freqs)
    prob_mat = np.zeros([n_vertex, n_freqs])

    # Loop across all vertices and subjects
    for vertex in range(n_vertex):
        for subj in range(n_subj):

            # Store centers for current vertex, current subj in temp vector
            cens_temp = centers[vertex, :, subj]

            # Loop through freq-ranges, counting when oscillations occur
            for ind, freq in enumerate(freqs):

                # Get the oscillation centers
                cens_fwin = _get_all_osc(cens_temp, freq, freq+f_win)

                # If there is an osc in range, add to prob_mat count
                if len(cens_fwin) != 0:
                    prob_mat[vertex, ind] += 1

    # Divide by # of subjects to get probability per freq-range
    prob_mat = prob_mat/n_subj

    # Initialize vectors to store correlations and p-values
    corr_vec = np.zeros([n_freqs-1])
    p_vec = np.zeros([n_freqs-1])

    # Compute corr between f and f+f_win start windows
    for f_ind in range(n_freqs-ind_step):
        corr_vec[f_ind], p_vec[f_ind] = pearsonr(prob_mat[:, f_ind], prob_mat[:, f_ind+ind_step])

    # Select frequency range that represents the start of each correlation
    freqs = freqs[:-1]

    return corr_vec, p_vec, freqs


def osc_space_group(oscs, bands, verts, osc_param=0, space_param=1):
    """

    Parameters
    ----------
    oscs : dict()
        xx
    bands : ?
        xx
    verts : 2d array
        Spatial coordinates for all vertices [n_verts, 3].
    osc_param : ?
        xx
    space_param : ?
        xx

    Returns
    -------
    dat_out : 3d array
        Correlation data for all subjects, all bands.
            [n_subjs, n_bands, 2], where last dimension is [R-val, p-val].
    labels : list of str
        Labels of oscillation bands that were analyzed in dat_out.
    """

    labels = list(bands.keys())
    n_verts, n_bands, n_subjs = oscs[labels[0]].shape

    space = verts[:, space_param]

    dat_out = np.zeros(shape=(n_subjs, len(bands), 2))

    for subj in range(n_subjs):
        for ind, band in enumerate(bands):

            cur_dat = oscs[band][:, osc_param, subj]
            keep_inds = cur_dat > 0

            dat_out[subj, ind, :] = pearsonr(cur_dat[keep_inds], space[keep_inds])

    return dat_out, labels

#################################################################################################
############################ OMEGAMAPPIN - OM_MD - PRIVATE FUNCTIONS ############################
#################################################################################################

def _get_all_osc(centers, osc_low, osc_high):
    """Returns all the oscillations in a specified frequency band.

    Parameters
    ----------
    centers : 1d array
        Vector of oscillation centers.
    osc_low : int
        Lower bound for frequency range.
    osc_high : int
        Upper bound for frequency range.

    Returns
    -------
    osc_cens : 1d array
        Osc centers in specified frequency band.
    """

    # Get inds of desired oscs and pull out from input data
    osc_inds = (centers >= osc_low) & (centers <= osc_high)
    osc_cens = centers[osc_inds]

    return osc_cens


def _osc_prob(osc_mat):
    """Takes a 3D matrix of oscillations across subjects, calculates probability of oscillation.

    Parameters
    ----------
    osc_mat : 3d array
        Oscillations for each subject, [n_vertex, n_dim, n_subj].

    Returns
    -------
    prob : 1d array
        Vector with probability of given oscillation at each vertex.
    """

    # Check how many vertices and subjects in group
    [n_vertex, n_dim, n_subj] = np.shape(osc_mat)

    # Initialize vector to store probabilities
    prob = np.zeros([n_vertex])

    # Loop through all vertices, calculating osc prob for each
    for i in range(n_vertex):
        prob[i] = (np.count_nonzero(osc_mat[i, 0, :]) / n_subj)

    return prob


def _osc_pow_ratio(osc_mat):
    """Calculate the power ratio of an oscillation.

    Power ratio is a score, bounded between [0, 1], reflecting power
    in a given frequency band, relative to the max power in that
    frequency band.

    Parameters
    ----------
    osc_mat : 3d array
        Oscillations for each subject, [n_vertex, n_dim, n_subj].

    Returns
    -------
    pow_ratio : 1d array
        Vector with oscillation score of given oscillation at each vertex.
    """

    # Check how many vertices and subjects in group
    [n_vertex, n_dim, n_subj] = np.shape(osc_mat)

    # Initialize vector to store average powers
    avg_powers = np.zeros(n_vertex)

    # Loop through all vertices
    for vertex in range(n_vertex):

        # Pull out temp vector of all oscillation powers
        temp_pows = osc_mat[vertex, 1, :]
        temp_pows = temp_pows[np.nonzero(temp_pows)]

        # If there are oscillations get average power
        if len(temp_pows) == 0:
            avg_powers[vertex] = 0
        else:
            avg_powers[vertex] = np.mean(temp_pows)

    # Get the maximum power across all vertices
    max_all = max(avg_powers)

    # Initialize vector to store power ratios
    pow_ratio = np.zeros(n_vertex)

    # Loop through all vertices, calculating power ratio
    for vertex in range(n_vertex):
        pow_ratio[vertex] = np.mean(osc_mat[vertex, 1, :]) / max_all

    return pow_ratio


def _band_sort(osc_bands):
    """Sort oscillation dictionary into order.

    Parameters
    ----------
    osc_bands : dict
        A dictionary containing the oscillation band definitions.

    Returns
    -------
    ordered_bands : list of str
        A list of the oscillation band names, in order.
    sort_inds : list of int
        A list of indices to sort oscillation bands.
    """

    # Check how many oscillation bands there are
    n_bands = len(osc_bands)

    # Initialize to store names and lower bounds
    band_names = []
    lower_bounds = np.array([])

    # Loop through, and grab name and lower bound for each band
    for band in osc_bands:
        band_names.append(band)
        lower_bounds = np.append(lower_bounds, osc_bands[band][0])

    # Get the indices to sort the lower bounds
    sort_inds = np.argsort(lower_bounds)

    # Use indices to sort band names into order
    ordered_bands = []
    ordered_bands[:] = [band_names[i] for i in sort_inds]

    return ordered_bands, sort_inds
