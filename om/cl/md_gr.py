"""MEG-DATA (MD) Analysis Module - Group

This code...

"""

from __future__ import print_function, division
import os
import pickle
import datetime
import numpy as np
import scipy.io as sio
from scipy.stats.stats import pearsonr

# Import custom om code
from om.gen import DataNotComputedError, InconsistentDataError, UnknownDataTypeError

##########################################################################################
############################  OMEGAMAPPIN - MD_GROUP CLASSES  ############################
##########################################################################################

class GroupMegData(MegData):
    """A class to store OMEGA data from multiple subjects.

    Holds all oscillations, regardless of spatial location.
    Note: Class derived from MegData()
    """

    def __init__(self, db, osc):
        """

        Parameters
        ----------
        db : ?
            xx
        osc : ?
            xx
        """

        # Initialize from MegData() object
        MegData.__init__(self, db)

        # Initialize groups subject variables
        self.n_subjs = int()
        self.subjs = []

        # Set definition of oscillation bands used for the group
        self.bands = osc.bands

        # Initialize count of total oscillations, across all subjects
        self.n_oscs_tot = int()

        # Set title for plots
        self.title = 'Group'

        # Initialize dictionary for oscillation band data - NEW
        self.gr_oscs = dict()

        # Initilaize dictionary to store oscillation probabilities - NEW
        self.osc_probs = dict()

        # Initialize dict to store oscillation power ratios - NEW
        self.osc_pow_ratios = dict()

        # Initialize to store oscillation scores - NEW
        self.osc_scores = dict()

        # Initialize vars to store slope values
        self.vert_slopes = np.array([])
        self.slopes_gr_avg = np.array([])

        # Set booleans for what has been run
        self.osc_prob_done = False
        self.osc_score_done = False


    def add_subject(self, new_subj, add_all_oscs=False, add_vertex_bands=False,
                    add_vertex_oscs=False, add_peak_freqs=False, add_vertex_slopes=False):
        """Adds a new subject to the GroupMegData object.

        Parameters
        ----------
        self : GroupMegData() object.
            Object to store map data across a group of subjects.
        new_subj : MegData() Object
            MEG subject (instance of MegData)
        add_all_oscs : boolean, optional (default: False)
            Whether to add the vectors of all oscillations, collapsed across vertices.
        add_vertex_bands : boolean, optional (default: False)
            Whether to add the oscillation band data, across vertices.
        add_vertex_oscs : boolean, optional (default: False)
            Whether to add all oscillations, across vertices.
        add_peak_freqs : boolean, optional (default: False)
            Whether to add peak frequencies.
        add_vertex_slopes : boolean, optional (default: False)
            Whether to add the slopes.
        """

        # Check if subject has data
        if not new_subj.has_data:
            print("Empty meg data object. Cannot add data.")

        # Add All-Osc Data
        if add_all_oscs:

            # Add oscillation parameters to current data
            self.centers_all = np.append(self.centers_all, new_subj.centers_all)
            self.bws_all = np.append(self.bws_all, new_subj.bws_all)
            self.powers_all = np.append(self.powers_all, new_subj.powers_all)
            self.slopes = np.append(self.slopes, new_subj.slopes)

            # Add centers hist
            self.centers_hist.append(new_subj.centers_hist)

            # Update count of total number of oscillations
            self.n_oscs = np.append(self.n_oscs, new_subj.n_oscs)
            self.n_oscs_tot = len(self.centers_all)

        # Add band-specific data
        if add_vertex_bands:

            # Check that new subject has same bands defined
            if not set(self.bands) == set(new_subj.bands):
                raise InconsistentDataError('Oscillation bands are inconsistent.')

            # Add new subject to group oscillations
            if self.n_subjs == 0:
                self.gr_oscs = new_subj.oscs
            else:
                for band in self.bands:
                    self.gr_oscs[band] = np.dstack([self.gr_oscs[band], new_subj.oscs[band]])

        # Add oscillations per vertex
        if add_vertex_oscs:

            if self.n_subjs == 0:
                self.centers = new_subj.centers
                self.powers = new_subj.powers
                self.bws = new_subj.bws
            else:
                self.centers = np.dstack([self.centers, new_subj.centers])
                self.powers = np.dstack([self.powers, new_subj.powers])
                self.bws = np.dstack([self.bws, new_subj.bws])

        # Add oscillation peak data
        if add_peak_freqs:

            # Check that new subject has same bands defined
            if not set(self.bands) == set(new_subj.bands):
                raise InconsistentDataError('Oscillation bands are inconsistent.')

            #
            if self.n_subjs == 0:
                self.peaks = new_subj.peaks
            else:
                for band in self.bands:
                    self.peaks[band] = np.append(self.peaks[band], new_subj.peaks[band])

        # Add slopes per vertex
        if add_vertex_slopes:

            if self.n_subjs == 0:
                self.vert_slopes = new_subj.slopes
            else:
                self.vert_slopes = np.hstack([self.vert_slopes, new_subj.slopes])

        # Update subj count and subject number list
        self.n_subjs += 1
        self.subjs = np.append(self.subjs, new_subj.subnum)

        # Update booleans about what is loaded
        self.all_osc = add_all_oscs
        self.bands_vertex = add_vertex_bands

        # Add demographic data
        self.sex.append(new_subj.sex)
        self.age = np.append(self.age, new_subj.age)


    def group_slope(self, avg='mean'):
        """Calculates the average slope value for each vertex, across subjects.

        Parameters
        ----------
        self : GroupMegData() object.
            Object to store map data across a group of subjects.
        avg : {'mean', 'median'}, optional
            How to average across the group.
        """

        # Calculate the average slope value per vertex
        if avg is 'mean':
            self.slopes_gr_avg = np.mean(self.vert_slopes, axis=1)
        elif avg is 'median':
            self.slopes_gr_avg = np.median(self.vert_slopes, axis=1)


    def osc_prob(self):
        """Calculates the probability of an osc in a specific band.

         This is done per vertex, across subjects.
         """

        # Check if vertex data is set
        if not self.bands_vertex:
            raise DataNotComputedError('Vertex oscillation bands data not available.')

        # For each oscillation band, compute the probability of an oscillation in that band - NEW
        for band in self.bands:
            self.osc_probs[band] = _osc_prob(self.gr_oscs[band])

        # Update boolean that oscillation probability has been computed
        self.osc_prob_done = True


    def osc_score(self):
        """Calculate the oscillation score for each frequency band.

        The oscillation score is ....
        """

        # Check if oscillation probability is calculated. Can't proceed if it isnt.
        if not self.osc_prob_done:
            raise DataNotComputedError('Oscillation probability not computed - can not proceed.')

        # Compute power ratio for each oscillation band - NEW
        for band in self.bands:
            self.osc_pow_ratios[band] = _osc_pow_ratio(self.gr_oscs[band])

        # Compute oscillation score for each oscillation band - NEW
        for band in self.bands:
            self.osc_scores[band] = self.osc_pow_ratios[band] * self.osc_probs[band]

        # Set boolean that oscillation score has been computed.
        self.osc_score_done = True


    def osc_map_corrs(self, map_type):
        """Calculates the correlations between oscillation probabilities or scores.

        Parameters
        ----------
        self : GroupMegData() object.
            Object to store map data across a group of subjects.
        map_type : {'prob', 'score'}
            Which map data type to save out.

        Returns
        -------
        corrs : dict
            Contains the correlations between all oscillation probabilities.
        corrs_mat : array
            xx
        """

        # Check if oscillation probabilities have been calculated.
        if not self.osc_prob_done:
            raise DataNotComputedError('Oscillation probability not computed - can not proceed.')

        # Check how many oscillation bands are defined
        n_bands = len(self.bands)

        # Initialize matrices to store correlation results
        corrs_mat = np.zeros([n_bands, n_bands])
        ps_mat = np.zeros([n_bands, n_bands])

        # Get oscillation bands in order
        sorted_bands, sort_inds = _band_sort(self.bands)

        # Set which map to run
        if map_type is 'prob':
            dat = self.osc_probs
        elif map_type is 'score':
            dat = self.osc_scores
        else:
            raise UnknownDataTypeError('Map type not understood.')

        # Loop through all bands, computing correlations between them
        for i in range(n_bands):
            for j in range(n_bands):
                corrs_mat[i, j], ps_mat[i, j] = pearsonr(
                    dat[sorted_bands[sort_inds[i]]],
                    dat[sorted_bands[sort_inds[j]]])

        # Set diagonals to zero - where band is correlated with itself
        np.fill_diagonal(corrs_mat, 0)
        np.fill_diagonal(ps_mat, 0)

        return corrs_mat, ps_mat, sorted_bands


    def calc_osc_peak_age(self):
        """Compares age and peak frequency within frequency bands.

        Returns
        -------
        corrs_mat : dict
            A dictionary containing correlations results comparing age to oscillations.
        ps_mat : ?
            xx
        sorted_bands : ?
            xx
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


    def freq_corr(self, f_win):
        """Calculates the correlation between adjacent frequency bands.

        Uses oscillation probabilities.

        Parameters
        ----------
        self : GroupMegData() object.
            Object to store map data across a group of subjects.
        f_win : float
            Size of frequency window to use.

        Returns
        -------
        corr_vec : 1d array
            Vector of the correlation coefficients between all adjacent frequency bands.
        p_vec : 1d array
            Vector of the p-values for the correlations between adjacent frequency bands.
        """

        # Get # vertices, # of subjects to loop through
        [n_vertex, n_slots, n_subj] = np.shape(self.centers)

        # Initialize variables for # of freqs, and matrix to store probability
        n_freqs = len(range(3, 40-f_win))
        prob_mat = np.zeros([n_vertex, n_freqs])

        # Loop through all vertices
        for vertex in range(0, n_vertex):

            # Loop through all subjects
            for subj in range(0, n_subj):

                # Store centers for current vertex, current subj in temp vector
                cens_temp = self.centers[vertex, :, subj]

                # Loop through freq-ranges, counting when oscillations occur
                i = 0
                for freq in range(3, 40-f_win):

                    # Get the oscillation centers
                    cens_fwin = _get_all_osc(cens_temp, freq, freq + f_win)

                    # If there is an osc in range, add to prob_mat count
                    if len(cens_fwin) != 0:
                        prob_mat[vertex, i] += 1

                    i += 1

        # Divide by # of subjects to get probability per freq-range
        prob_mat = prob_mat/n_subj

        # Initialize vectors to store correlations and p-values
        corr_vec = np.zeros([n_freqs-1])
        p_vec = np.zeros([n_freqs-1])

        # Compute corr between f and f+f_win start windows
        for f_ind in range(0, n_freqs-f_win):
            corr_vec[f_ind], p_vec[f_ind] = pearsonr(prob_mat[:, f_ind], prob_mat[:, f_ind+f_win])

        return corr_vec, p_vec


    def save_gr_slope(self, file_name):
        """Saves out the average group slope results.

        Parameters
        ----------
        self : GroupMegData() object.
            Object to store map data across a group of subjects.
        file_name : str
            File name to save group slope file as.
        """

        # Set up
        pickle_file_name = file_name + '.p'
        pickle_save_name = os.path.join(self.db.maps_path, 'Slopes', pickle_file_name)

        # Check current time for when file is saved
        cur_time = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # Collect data together to save out
        dat_out = dict({'dat_source': self.dat_source,
                        'slopes': self.slopes_gr_avg,
                        'n_subjs': self.n_subjs,
                        'save_time': cur_time})

        # Save out with pickle
        pickle.dump(dat_out, open(pickle_save_name, 'wb'))


    def save_map(self, map_type, file_name):
        """Save oscillation map data out to disc.

        Parameters
        ----------
        self : GroupMegData() object.
            Object to store map data across a group of subjects.
        map_type : {'prob', 'score'}
            Which map data type to save out.
        file_name : str
            String to add to the file name.
        """

        # Set data type
        if map_type is 'prob':
            dat = self.osc_probs
        elif map_type is 'score':
            dat = self.osc_scores
        else:
            raise UnknownDataTypeError('Map type not understood.')

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


    def set_slope_viz(self):
        """Saves out a matfile, of the group average slope, for visualization."""

        # Set up paths to save to
        save_name = 'Group_Slopes'
        save_file = os.path.join(self.db.viz_path, save_name)

        # Save desired outputs into a dictionary
        save_dict = {}
        save_dict['slopes'] = self.slopes_gr_avg
        save_dict['dat_source'] = self.dat_source
        save_dict['n_subjs'] = self.n_subjs
        save_dict['save_time'] = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # Save the dicionary out to a .mat file
        sio.savemat(save_file, save_dict)


    def set_map_viz(self, map_type):
        """Set an oscillation map for visualization with Brainstorm.

        Parameters
        ----------
        map_type : {'prob', 'score'}
            Which map data type to set as viz.
        """

        # Set data type
        if map_type is 'prob':
            save_name = 'Group_Osc_Prob_Viz'
            dat = self.osc_probs
        elif map_type is 'score':
            dat = self.osc_scores
            save_name = 'Group_Osc_Score_Viz'
        else:
            raise UnknownDataTypeError('Map type not understood.')

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


#################################################################################################
############################ OMEGAMAPPIN - OM_MD - PRIVATE FUNCTIONS ############################
#################################################################################################

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
    for i in range(0, n_vertex):
        prob[i] = (np.count_nonzero(osc_mat[i, 0, :]) / n_subj)

    return prob


def _osc_pow_ratio(osc_mat):
    """Calculate the power ratio of an oscillation.

    Power ratio is a score, bounded between [0, 1], reflecting power
    in a given freqeuncy band, relative to the max power in that
    frequency band. Max power is ...

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
    for vertex in range(0, n_vertex):

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
    for vertex in range(0, n_vertex):
        pow_ratio[vertex] = np.mean(osc_mat[vertex, 1, :]) / max_all

    return pow_ratio


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
    osc_inds = (centers > osc_low) & (centers < osc_high)
    osc_cens = centers[osc_inds]

    return osc_cens


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
