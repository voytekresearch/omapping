from __future__ import print_function, division
import os
import csv
import pickle
import datetime
import numpy as np
import scipy.io as sio
from scipy.stats.stats import pearsonr
from om.gen import *

##########################################################################################
###############################  OMEGAMAPPIN - MD CLASSES  ###############################
##########################################################################################


class MegData(object):
    """Class for a single subject of FOOF results for MEG Source PSDs."""

    def __init__(self, db):
        """

        Parameters
        ----------
        db : ?
            xx
        """

        # Store which db is set
        self.dat_source = db.dat_source

        # Add database object
        self.db = db

        # Pull out needed paths from OMDB object
        #self.project_path = db.project_path
        #self.maps_path = db.maps_path
        #self.meg_path = db.meg_path
        #self.foof_path = db.foof_path
        #self.viz_path = db.viz_path
        #self.md_save_path = db.md_save_path

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

        # Initialize dictionary to store osc bands - NEW
        self.bands = dict()

        # Initialize dictionary for oscillation bands - NEW
        self.oscs = dict()

        # Initialize dictionary for peak frequencies - NEW
        self.peaks = dict()

        # Set plot title
        self.title = ''

        # Set boolean for what has been run
        self.has_data = False
        self.all_osc = False
        self.bands_vertex = False

        # Initialize demographic variables
        self.sex = []
        self.age = np.array([])

        # Initialize oscillation count
        self.osc_count = int()


    def import_foof(self, subnum, get_demo=True, load_type='pickle'):
        """Import FOOF results to MegData object.

        Parameters
        ----------
        self : MegData() object
            Object to store MEG FOOF data.
        subnum : int
            Number of the subject to import.
        get_demo : boolean, optional (default = True)
            Whether to load demographic data from csv file.
        load_type : {'pickle', 'csv'}, optional
            What type of file to load data from.
        """

        # Check if object already has data
        if self.has_data:
            print("Subject object already contains data. Can't add")
            return

        # Set subject number for current data object
        self.subnum = subnum
        self.title = 'S-' + str(self.subnum)

        # Set up paths, get list of files for available subjects
        files = os.listdir(os.path.join(self.db.foof_path, load_type))
        files = clean_file_list(files, 'Foof_Vertex')

        # Get specific file path for specific subject
        cur_subj_file = get_cur_subj(subnum, files)
        cur_subj_path = os.path.join(self.db.foof_path, load_type, cur_subj_file)

        # Load data file
        if load_type is 'pickle':
            results = _load_foof_pickle(cur_subj_path)
        elif load_type is 'csv': # NOTE: not yet implemented
            results = _load_foof_csv(cur_subj_path)

        # Pull out data from results - NOTE: New version.
        self.centers, self.powers, self.bws, self.slopes, self.n_psds \
            = extract_foof_pickle(results)

        """OLD:
        # Check how many psds there are
        self.n_psds = len(results)

        # Initialize numpy arrays to pull out different result params
        self.slopes = np.zeros([self.n_psds, 1])
        self.centers = np.zeros([self.n_psds, 8])
        self.powers = np.zeros([self.n_psds, 8])
        self.bws = np.zeros([self.n_psds, 8])

        # Loop through FOOF results, pulling out individual findings
        for i in range(self.n_psds):
            self.slopes[i] = results[i][0]
            self.centers[i, 0:len(results[i][1])] = results[i][1]
            self.powers[i, 0:len(results[i][2])] = results[i][2]
            self.bws[i, 0:len(results[i][3])] = results[i][3]
        """

        # Check how many oscillations per vertex
        self.osc_count = np.zeros([self.n_psds, 1])
        for i in range(0, self.n_psds):
            self.osc_count[i, 0] = len(np.nonzero(self.centers[i, :])[0])

        # Get demographic data
        if get_demo:
            self.sex, self.age = _get_demo_csv(self.subnum, self.db.meg_path, self.dat_source)

        # Update boolean to say current subject has data attached
        self.has_data = True


    def osc_bands_vertex(self, osc, avg='mean'):
        """Groups oscillations at each vertex in distinct frequency bands.
        Stores band specific oscillations in (self.){thetas, alphas, betas, lowgammas}.

        Parameters
        ----------
        self : MegData object
            MegData object.
        osc : Osc object
            An object containing frequency bands to use.
        avg : {'mean', 'median'}, optional
            How to average to calculate peak frequencies.
        """

        # Save bands used
        self.bands = osc.bands

        # Initialize matrices to store oscillations in each band
        for band in osc.bands:
            self.oscs[band] = np.zeros([self.n_psds, 4])

        # Loop through each vertex
        for i in range(self.n_psds):

            # Get centers, powers and bws from individual vertex
            centers_temp = self.centers[i, :]
            powers_temp = self.powers[i, :]
            bws_temp = self.bws[i, :]

            # Loop through each oscillation band
            for band in osc.bands:

                # Get oscillations in specific band
                self.oscs[band][i, :] = _get_single_osc(centers_temp, powers_temp, bws_temp,
                                                        osc.bands[band][0], osc.bands[band][1])

        # Update boolean to note that current subject has band specific oscs calculated.
        self.bands_vertex = True


    def all_oscs(self, verbose=True):
        """Flatten osc data to vectors.

        Parameters
        ----------
        verbose : xx
            xx

        When imported, oscillation data is in matrix form [n_vertex, osc_slots].
        This functions converts these matrices into 1-D vectors.

        Note: This function loses information about which vertex oscillations occur at.
        """

        # Flatten osc data into vectors. Uses C-style row-major order
        self.centers_all = self.centers.flatten('C')
        self.powers_all = self.powers.flatten('C')
        self.bws_all = self.bws.flatten('C')

        # Flattened vectors will have lots of zeros. Get only non-zero indices.
        non_zeros = np.nonzero(self.centers_all)
        self.centers_all = self.centers_all[non_zeros]
        self.powers_all = self.powers_all[non_zeros]
        self.bws_all = self.bws_all[non_zeros]

        # Check for nans in BW estimation
        # NOTE: Updated FOOF sometimes returns NaN for bw. Check and discard those.
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
        self.all_osc = True


    def peak_freq(self, osc, avg='mean'):
        """Calculates the peak frequency for each oscillatory band.

        Parameters
        ----------
        osc : Osc object
            Object with oscillation frequency details
        avg : {'mean', 'median'}, optional
            Which type of averaging to do.
        """

        # Loop through each band, calculating peak frequency
        for band in osc.bands:
            self.peaks[band] = _osc_peak(
                self.centers_all, osc.bands[band][0], osc.bands[band][1], avg)


    def calc_osc_param_corrs(self):
        """Calculates correlations between oscillatory parameters."""

        # Set labels for the things being correlated
        labels = ['Centers', 'Powers', 'Bandwidths']

        # Check how many categories there are
        n = len(labels)

        # Initialize matrices to store R and p values
        corrs_mat = np.zeros([n, n])
        ps_mat = np.zeros([n, n])

        # Calculate correlations between all parameters
        corrs_mat[0, 1], ps_mat[0, 1] = pearsonr(self.centers_all, np.log10(self.powers_all))
        corrs_mat[0, 2], ps_mat[0, 2] = pearsonr(self.centers_all, np.log10(self.bws_all))
        corrs_mat[1, 2], ps_mat[1, 2] = pearsonr(np.log10(self.powers_all),
                                                 np.log10(self.bws_all))

        # Duplicate matrix across the diagonal
        corrs_mat = corrs_mat + corrs_mat.T
        ps_mat = ps_mat + ps_mat.T

        return corrs_mat, ps_mat, labels


    def set_foof_viz(self):
        """Saves a matfile of freq info to be loaded with Brainstorm for visualization."""

        # Set up paths to save to
        save_name = str(self.subnum) + '_Foof_Viz'
        save_file = os.path.join(self.db.viz_path, save_name)

        # Initialzie dictionary, save basic information and slope data
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


############################################################################################
########################## OMEGAMAPPIN - OM_MD - PUBLIC FUNCTIONS ##########################
############################################################################################

def print_corrs_mat(rs_mat, ps_mat, labels):
    """Prints out correlations - from a given matrix.

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
    n = len(labels)

    # Loop through the matrix to print out
    for x in range(n):
        for y in range(n):

            # Skip bottom triangle and diagonal
            if x == y or y < x:
                continue

            # Print out correlation
            print('Corr of ', '{:18}'.format(labels[x]+'-'+labels[y]),
                  ' is ', '{:+1.4f}'.format(rs_mat[x, y]), '    with p-val of ',
                  '{:1.5f}'.format(ps_mat[x, y]))


def print_corrs_vec(rs_vec, ps_vec, labels, desc):
    """Prints out corrlations - from a given vector.

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
    n = len(labels)

    # Loop through vectors, printing out the correlations.
    for x in range(n):
        print('Corr of ', '{:20}'.format(labels[x]+'-'+desc), ' is ',
              '{:+1.4f}'.format(rs_vec[x]), '    with p-val of ',
              '{:1.5f}'.format(ps_vec[x]))


def save_md_pickle(obj, save_name):
    """Save current meg data object as a pickled object.

    Parameters
    ----------
    obj : MegData() or GroupMegData()
        Object to save to pickle
    save_name : str
        String to be included in the name of the file.
    """

    # Get database object
    db = OMDB()

    # Set save name and path
    save_name = 'Res_' + save_name + '_' + datetime.datetime.now().strftime("%Y-%m-%d") + '.p'
    foof_save_path = os.path.join(db.md_save_path, save_name)

    # Save out data to pickle file
    pickle.dump(obj, open(foof_save_path, 'wb'))


def load_md_pickle(file_name):
    """Load a pickled file.

    Parameters
    ----------
    file_name : str
        File name of the pickle file to be loaded.
    """

    # Get database object
    db = OMDB()

    # Set up path to load file
    dat_path = os.path.join(db.md_save_path, file_name)

    # Load file & return pickled object
    results = pickle.load(open(dat_path, 'rb'))

    return results


#################################################################################################
############################ OMEGAMAPPIN - OM_MD - PRIVATE FUNCTIONS ############################
#################################################################################################


def _get_single_osc(centers, powers, bws, osc_low, osc_high):
    """ Searches for an oscillations of specified frequency band.

    Returns a single oscillation in that band.
    Helper function for osc_per_vertex in MegData.

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
    osc_out : tuple
        Osc data, form: [centers, powers, bws, # oscillations].
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


def _get_demo_csv(subnum, meg_path, dat_source):
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
        csv_file_name = '00-OMEGA_Subjects.csv'
        num_ind = 1
        sex_ind = 4
        age_ind = 7
    elif dat_source is 'HCP':
        csv_file_name = '00-HCP_Subjects.csv'
        num_ind = 0
        sex_ind = 3
        age_ind = 4
    csv_file = os.path.join(meg_path, csv_file_name)

    # Open csv file, loop through looking for right row, grab age & sex information
    with open(csv_file, 'rb') as f_name:
        reader = csv.reader(f_name, delimiter=',')
        for row in reader:
            if row[num_ind] == str(subnum):
                sex = row[sex_ind]
                if dat_source is 'OMEGA':
                    age = int(row[age_ind])
                else:
                    age_temp = (row[age_ind]).split('-')
                    age = (int(age_temp[0]) + int(age_temp[1]))/2
                break

    return sex, age


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


def _osc_peak(centers, osc_low, osc_high, avg='mean'):
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
        Peak frequency value - the average frequency within a given range.
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


def _load_foof_pickle(file_name):
    """Loads FOOF data from a pickle file.

    Parameters
    ----------
    filename : str
        Full path, including filename, to file to be loaded.
    """

    # Load from pickle file
    results = pickle.load(open(file_name, 'rb'))

    return results


def _load_foof_csv(file_name):
    """
    NOTE: not yet implemented
    """

    pass
