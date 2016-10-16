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


class MegData():
    """Class for a single subject of FOOF results for MEG Source PSDs."""

    def __init__(self, db):

        # Store which db is set
        self.dat_source = db.dat_source

        # Pull out needed paths from OMDB object
        self.project_path = db.project_path
        self.maps_path = db.maps_path
        self.meg_path = db.meg_path
        self.foof_path = db.foof_path
        self.viz_path = db.viz_path
        self.md_save_path = db.md_save_path

        # Initialize subject number
        self.subnum = int()
        self.n_PSDs = int()

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


        # Initialize arrays for oscillation bands - OLD
        self.thetas = np.array([])
        self.alphas = np.array([])
        self.betas = np.array([])
        self.lowgammas = np.array([])

        # Initialize peak frequency variables - OLD
        self.peak_theta = np.array([])
        self.peak_alpha = np.array([])
        self.peak_beta = np.array([])
        self.peak_lowgamma = np.array([])


    def import_foof(self, subnum, get_demo=True, load_type='pickle'):
        """Import FOOF results to MegData object.

        Parameters
        ----------
        self : MegData() object
            Object to store MEG FOOF data.
        subnum : int
            Number of the subject to import.
        get_demo : boolean, optional
            Whether to load demographic data from csv file.
        load_type : str, optional
            What type of file to load data from.
                Options: 'pickle', 'csv'
        """

        # Check if object already has data
        if self.has_data:
            print("Subject object already contains data. Can't add")
            return

        # Set subject number for current data object
        self.subnum = subnum
        self.title = 'S-' + str(self.subnum)

        # Set up paths, get list of files for available subjects
        files = os.listdir(os.path.join(self.foof_path, load_type))
        files = clean_file_list(files, 'Foof_Vertex')

        # Get specific file path for specific subject
        cur_subj_file = get_cur_subj(subnum, files)
        cur_subj_path = os.path.join(self.foof_path, load_type, cur_subj_file)

        # Load data file
        if load_type is 'pickle':
            results = _load_foof_pickle(cur_subj_path)
        elif load_type is 'csv':
            results = _load_foof_csv(cur_subj_path)

        # Check how many
        self.n_PSDs = len(results)

        # Initialize numpy arrays to pull out different result params
        self.slopes = np.zeros([self.n_PSDs, 1])
        self.centers = np.zeros([self.n_PSDs, 8])
        self.powers = np.zeros([self.n_PSDs, 8])
        self.bws = np.zeros([self.n_PSDs, 8])

        # Loop through FOOF results, pulling out individual findings
        for i in range(0, self.n_PSDs):
            self.slopes[i] = results[i][0]
            self.centers[i, 0:len(results[i][1])] = results[i][1]
            self.powers[i, 0:len(results[i][2])] = results[i][2]
            self.bws[i, 0:len(results[i][3])] = results[i][3]

        # Check how many oscillations per vertex
        self.osc_count = np.zeros([self.n_PSDs, 1])
        for i in range(0, self.n_PSDs):
            self.osc_count[i, 0] = len(np.nonzero(self.centers[i, :])[0])

        # Get demographic data
        if get_demo:
            self.sex, self.age = _get_demo_csv(self.subnum, self.meg_path, self.dat_source)

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
        avg : str, optional
            How to average to calculate peak frequencies.
                Options: {'mean', 'median'}
        """
        """
        ## Re-Initialize matrices to right size to save results
        self.thetas = np.zeros([self.n_PSDs, 4])
        self.alphas = np.zeros([self.n_PSDs, 4])
        self.betas = np.zeros([self.n_PSDs, 4])
        self.lowgammas = np.zeros([self.n_PSDs, 4])

        # Loop through each vertex
        for i in range(0, self.n_PSDs):

            # Get centers, powers and bws from individual vertex
            centers_temp = self.centers[i, :]
            powers_temp = self.powers[i, :]
            bws_temp = self.bws[i, :]

            # Get oscillations in specific band for each band
            self.thetas[i, :]    = _get_osc(centers_temp, powers_temp, bws_temp,
                                           osc.theta_low, osc.theta_high)
            self.alphas[i, :]    = _get_osc(centers_temp, powers_temp, bws_temp,
                                           osc.alpha_low, osc.alpha_high)
            self.betas[i, :]     = _get_osc(centers_temp, powers_temp, bws_temp,
                                           osc.beta_low, osc.beta_high)
            self.lowgammas[i, :] = _get_osc(centers_temp, powers_temp, bws_temp,
                                           osc.lowgamma_low, osc.lowgamma_high)

        # Update boolean to note that current subject has band specific oscs calculated.
        self.bands_vertex = True
        """

        ## NEW VERSION

        # Save bands used
        self.bands = osc.bands

        # Initialize matrices to store oscillations in each band
        for band in osc.bands:
            self.oscs[band] = np.zeros([self.n_PSDs, 4])

        # Loop through each vertex
        for i in range(self.n_PSDs):

            # Get centers, powers and bws from individual vertex
            centers_temp = self.centers[i, :]
            powers_temp = self.powers[i, :]
            bws_temp = self.bws[i, :]

            # Loop through each oscillation band
            for band in osc.bands:

                # Get oscillations in specific band
                self.oscs[band][i, :] = _get_osc(centers_temp, powers_temp, bws_temp,
                                                 osc.bands[band][0], osc.bands[band][1])

        # Update boolean to note that current subject has band specific oscs calculated.
        self.bands_vertex = True


    def save_viz(self):
        """Saves a matfile of freq info to be loaded with Brainstorm for visualization.
        NOTE: Needs updating for Osc_Dict
        """

        # Set up paths to save to
        save_name = str(self.subnum) + '_Foof_Viz'
        save_file = os.path.join(self.viz_path, save_name)

        # Save desired outputs into a dictionary
        save_dict = {}
        save_dict['slopes'] = self.slopes
        save_dict['thetas'] = self.thetas
        save_dict['alphas'] = self.alphas
        save_dict['betas'] = self.betas
        save_dict['lowgammas'] = self.lowgammas

        # Save the dicionary out to a .mat file
        sio.savemat(save_file, save_dict)


    def all_oscs(self, verbose=True):
        """Flatten osc data to vectors.

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
        avg : str
            Which type of averaging to do.
                Options: {'mean', 'median'}
        """

        """
        # Get peak frequency within each frequency band
        self.peak_theta = _osc_peak(self.centers_all, osc.theta_low, osc.theta_high, avg)
        self.peak_alpha = _osc_peak(self.centers_all, osc.alpha_low, osc.alpha_high, avg)
        self.peak_beta = _osc_peak(self.centers_all, osc.beta_low, osc.beta_high, avg)
        self.peak_lowgamma = _osc_peak(self.centers_all, osc.lowgamma_low, osc.lowgamma_high, avg)
        """

        ## NEW VERSION

        # Loop through each band, calculating peak frequency
        for band in osc.bands:
            self.peaks[band] = _osc_peak(self.centers_all, osc.bands[band][0], osc.bands[band][1], avg)


    def calc_osc_param_corrs(self):
        """Calculates correlations between oscillatory parameters."""

        """
        # Calculate correlations: Centers vs. Bandwidth
        corr_cen_bw, pcorr_cen_bw = pearsonr(self.centers_all, np.log10(self.bws_all))

        # Calculate correlations: Centers vs. Power
        corr_cen_pow, pcorr_cen_pow = pearsonr(self.centers_all, np.log10(self.powers_all))

        # Calculate correlations: Bandwidth vs. Power
        corr_bw_pow, pcorr_bw_pow = pearsonr(np.log10(self.bws_all), np.log10(self.powers_all))

        # Save correlation results to list to return
        corrs = [['Center - B.W.'   , corr_cen_bw , pcorr_cen_bw ],
                 ['Center - Power'  , corr_cen_pow, pcorr_cen_pow],
                 ['B.W. - Power '   , corr_bw_pow , pcorr_bw_pow ]]

        return corrs
        """

        #### NEW

        labels = ['Centers', 'Powers', 'Bandwidths']

        n = len(labels)

        corrs_mat = np.zeros([n, n])
        ps_mat = np.zeros([n, n])

        #
        corrs_mat[0, 1], ps_mat[0, 1] = pearsonr(self.centers_all, np.log10(self.powers_all))
        corrs_mat[0, 2], ps_mat[0, 2] = pearsonr(self.centers_all, np.log10(self.bws_all))
        corrs_mat[1, 2], ps_mat[1, 2] = pearsonr(np.log10(self.powers_all),
                                                 np.log10(self.bws_all))

        corrs_mat = corrs_mat + corrs_mat.T
        ps_mat = ps_mat + ps_mat.T

        return corrs_mat, ps_mat, labels


class GroupMegData(MegData):
    """A class to store OMEGA data from multiple subjects.

    Holds all oscillations, regardless of spatial location.
    Note: Class derived from MegData
    """

    def __init__(self, db, osc):

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

        # Initialize matrices for oscillation band data - OLD
        self.gr_thetas = np.array([])
        self.gr_alphas = np.array([])
        self.gr_betas = np.array([])
        self.gr_lowgammas = np.array([])

        # Initialize dictionary for oscillation band data - NEW
        self.gr_oscs = dict()

        # Initialize to store oscillation probabilities - OLD
        self.theta_prob = np.array([])
        self.alpha_prob = np.array([])
        self.beta_prob = np.array([])
        self.lowgamma_prob = np.array([])

        # Initilaize dictionary to store oscillation probabilities - NEW
        self.osc_probs = dict()

        # Initialize to store oscillation power ratios - OLD
        self.theta_pow_ratio = np.array([])
        self.alpha_pow_ratio = np.array([])
        self.beta_pow_ratio = np.array([])
        self.lowgamma_pow_ratio = np.array([])

        # Initialize dict to store oscillation power ratios
        self.osc_pow_ratios = dict()

        # Initialize to store oscillation scores - OLD
        self.theta_score = np.array([])
        self.alpha_score = np.array([])
        self.beta_score = np.array([])
        self.lowgamma_score = np.array([])

        # Initialize to store oscillation scores - NEW
        self.osc_scores = dict()

        # Initialize vars to store slope values
        self.vert_slopes = np.array([])
        self.slopes_avg = np.array([])

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
        add_all_oscs : boolean, optional
            Whether to add the vectors of all oscillations, collapsed across vertices.
        add_vertex_bands : boolean, optional
            Whether to add the oscillation band data, across vertices.
        add_vertex_oscs : boolean, optional
            Whether to add all oscillations, across vertices.
        add_peak_freqs : boolean, optional
            Whether to add peak frequencies.
        add_vertex_slopes : boolean, optional
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

        """
        # Add band-specific data - OLD
        if add_vertex_bands:

            if self.n_subjs == 0:
                self.gr_thetas = new_subj.thetas
                self.gr_alphas = new_subj.alphas
                self.gr_betas = new_subj.betas
                self.gr_lowgammas = new_subj.lowgammas
            else:
                self.gr_thetas = np.dstack([self.gr_thetas, new_subj.thetas])
                self.gr_alphas = np.dstack([self.gr_alphas, new_subj.alphas])
                self.gr_betas = np.dstack([self.gr_betas, new_subj.betas])
                self.gr_lowgammas = np.dstack([self.gr_lowgammas, new_subj.lowgammas])
        """

        # Add band-specific data - NEW
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

        #
        if add_vertex_oscs:

            if self.n_subjs == 0:
                self.centers = new_subj.centers
                self.powers = new_subj.powers
                self.bws = new_subj.bws
            else:
                self.centers = np.dstack([self.centers, new_subj.centers])
                self.powers = np.dstack([self.powers, new_subj.powers])
                self.bws = np.dstack([self.bws, new_subj.bws])

        # Add oscillation peak data - NEW
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

        #
        if add_vertex_slopes:

            if self.n_subjs == 0:
                self.vert_slopes = new_subj.slopes
            else:
                self.vert_slopes = np.hstack([self.vert_slopes, new_subj.slopes])

        # Update subj count and subject number list
        self.n_subjs += 1
        self.subjs = np.append(self.subjs, new_subj.subnum)

        #
        self.all_osc = add_all_oscs
        self.bands_vertex = add_vertex_bands

        # Add demographic data
        self.sex.append(new_subj.sex)
        self.age = np.append(self.age, new_subj.age)

        """
        # Add oscillation peak data - OLD
        self.peak_theta = np.append(self.peak_theta, new_subj.peak_theta)
        self.peak_alpha = np.append(self.peak_alpha, new_subj.peak_alpha)
        self.peak_beta = np.append(self.peak_beta, new_subj.peak_beta)
        self.peak_lowgamma = np.append(self.peak_lowgamma, new_subj.peak_lowgamma)
        """


    def osc_prob(self):
        """Calculates the probability (per vertex / across subjects) of an osc in a specific band."""

        """"
        # For each frequency band, compute the probability of oscillation in that band. - OLD
        self.theta_prob = _osc_prob(self.gr_thetas)
        self.alpha_prob = _osc_prob(self.gr_alphas)
        self.beta_prob = _osc_prob(self.gr_betas)
        self.lowgamma_prob = _osc_prob(self.gr_lowgammas)
        """

        # Check if vertex data is set
        if not self.bands_vertex:
            raise DataNotComputedError('Vertex oscillation bands data not available.')

        # For each oscillation band, compute the probability of an oscillation in that band - NEW
        for band in self.bands:
            self.osc_probs[band] = _osc_prob(self.gr_oscs[band])

        # Update boolean that oscillation probability has been computed
        self.osc_prob_done = True


    def group_slope(self, save_out=False, file_name=None, set_viz=False):
        """Calculates the average slope value for each vertex, across subjects.

        NOTE: UPDATE CHIS SAVE NAME, NEED TO CHECK WHAT LOADS THESE FILES.

        Parameters
        ----------
        self : GroupMegData() object.
            Object to store map data across a group of subjects.
        save_out : boolean, optional
            Whether to save out npz file of average slope values.
        file_name : str
            File name to save group slope file as.
        set_viz : boolean, optional
            Whether to save out mat file for visualization in matlab.
        """

        # Calculate the average slope value per vertex
        self.slopes_avg = np.median(self.vert_slopes, axis=1)

        # Save map out, if required
        if save_out:

            # Set up
            npz_file_name = file_name + '.npz'
            npz_save_name = os.path.join(self.maps_path, 'Slopes', npz_file_name)

            # Save out an npz file
            np.savez(npz_save_name, chis=self.slopes_avg, n_subjs=self.n_subjs)

        # Save out as matfile for visualization with matlab, if required
        if set_viz:

            # Set up paths to save to
            save_name = 'Group_Slopes'
            save_file = os.path.join(self.viz_path, save_name)

            # Save desired outputs into a dictionary
            save_dict = {}
            save_dict['chis'] = self.slopes_avg
            save_dict['dat_source'] = self.dat_source
            save_dict['n_subjs'] = self.n_subjs
            save_dict['save_time'] = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

            # Save the dicionary out to a .mat file
            sio.savemat(save_file, save_dict)


    def set_prob_viz(self):
        """Saves out a matfile (of osc probs) to be loaded with Brainstorm for visualization. """

        # Set up paths to save to
        save_name = 'Group_Osc_Prob_Viz'
        save_file = os.path.join(self.viz_path, save_name)

        # Save desired outputs into a dictionary
        save_dict = {}
        save_dict['theta_prob'] = self.theta_prob
        save_dict['alpha_prob'] = self.alpha_prob
        save_dict['beta_prob'] = self.beta_prob
        save_dict['lowgamma_prob'] = self.lowgamma_prob
        save_dict['dat_source'] = self.dat_source
        save_dict['n_subjs'] = self.n_subjs
        save_dict['save_time'] = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # Save the dicionary out to a .mat file
        sio.savemat(save_file, save_dict)


    def osc_map_corrs(self, map_type):
        """Calculates the correlations between oscillation probabilities or scores.

        Returns
        -------
        corrs : dict
            Contains the correlations between all oscillation probabilities.
        corrs_mat : array
            xx
        """

        ## NEW

        # Check if oscillation probabilities have been calculated.
        if not self.osc_prob_done:
            raise DataNotComputedError('Oscillation probability not computed - can not proceed.')

        # Check how many oscillation bands are defined
        n_bands = len(self.bands)

        # Initialize matrices to store correlation results
        corrs_mat = np.zeros([n_bands, n_bands])
        ps_mat = np.zeros([n_bands, n_bands])

        # xx
        sorted_bands, sort_inds = _band_sort(self.bands)

        # xx
        if map_type is 'prob':
            dat = self.osc_probs
        elif map_type is 'score':
            dat = self.osc_scores
        else:
            raise UnknownDataTypeError('Map type not understood.')

        # xx
        for i in range(n_bands):
            for j in range(n_bands):
                corrs_mat[i, j], ps_mat[i, j] = pearsonr(
                    dat[sorted_bands[sort_inds[i]]],
                    dat[sorted_bands[sort_inds[j]]])

        # xx
        np.fill_diagonal(corrs_mat, 0)
        np.fill_diagonal(ps_mat, 0)

        return corrs_mat, ps_mat, sorted_bands

        """
        # Theta-Alpha Corr
        [r_th_al, p_th_al] = pearsonr(self.theta_prob, self.alpha_prob)

        # Theta-Beta Corr
        [r_th_be, p_th_be] = pearsonr(self.theta_prob, self.beta_prob)

        # Theta-LG Corr
        [r_th_lg, p_th_lg] = pearsonr(self.theta_prob, self.lowgamma_prob)

        # Alpha-Beta Corr
        [r_al_be, p_al_be] = pearsonr(self.alpha_prob, self.beta_prob)

        # Alpha-LG Corr
        [r_al_lg, p_al_lg] = pearsonr(self.alpha_prob, self.lowgamma_prob)

        # Beta-LG Corr
        [r_be_lg, p_be_lg] = pearsonr(self.beta_prob, self.lowgamma_prob)

        # Save correlation results in a list to return
        corrs = [['Theta-Alpha', r_th_al, p_th_al],
                 ['Theta-Beta' , r_th_be, p_th_be],
                 ['Theta-LG'   , r_th_lg, p_th_lg],
                 ['Alpha-Beta' , r_al_be, p_al_be],
                 ['Alpha-LG'   , r_al_lg, p_al_lg],
                 ['Beta-LG'    , r_be_lg, p_be_lg]]

        # Save corrs out in a matrix
        corrs_mat = np.array([[0, r_th_al, r_th_be, r_th_lg], [r_th_al, 0, r_al_be, r_al_lg],
                              [r_th_be, r_al_be, 0, r_be_lg], [r_th_lg, r_al_lg, r_be_lg, 0]])

        return corrs, corrs_mat
        """


    def osc_score(self):
        """Calculate the oscillation score for each frequency band.

        The oscillation score is .... XXXXX
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

        """
        # Compute power ratio for each oscillatory band - OLD
        self.theta_pow_ratio = _osc_pow_ratio(self.gr_thetas)
        self.alpha_pow_ratio = _osc_pow_ratio(self.gr_alphas)
        self.beta_pow_ratio = _osc_pow_ratio(self.gr_betas)
        self.lowgamma_pow_ratio = _osc_pow_ratio(self.gr_lowgammas)

        # Compute oscillation score for each oscillatory band - OLD
        self.theta_score = self.theta_pow_ratio * self.theta_prob
        self.alpha_score = self.alpha_pow_ratio * self.alpha_prob
        self.beta_score = self.beta_pow_ratio * self.beta_prob
        self.lowgamma_score = self.lowgamma_pow_ratio * self.lowgamma_prob
        """

        # Set boolean that oscillation score has been computed.
        self.osc_score_done = True


    def osc_score_corrs(self):
        """Calculates the correlations between oscillations scores.
        NOTE: NO LONGER NEEDED. NOW INCLUDED IN OSC_MAP_CORRS()

        Returns
        -------
        corrs : dict
            Contains the results of the correlations across oscillation scores.
        corrs_mat : array
            xx
        """

        ## NEW

        # Check if oscillation probabilities have been calculated.
        if not self.osc_score_done:
            raise DataNotComputedError('Oscillation score not computed - can not proceed.')

        # Check how many bands there are
        n_bands = len(self.bands)

        # Initialize matrices to store correlation results
        corrs_mat = np.zeros([n_bands, n_bands])
        ps_mat = np.zeros([n_bands, n_bands])

        # xx
        sorted_bands, sort_inds = _band_sort(self.bands)

        # xx
        for i in range(n_bands):
            for j in range(n_bands):
                corrs_mat[i, j], ps_mat[i, j] = pearsonr(
                    self.osc_scores[sorted_bands[sort_inds[i]]],
                    self.osc_scores[sorted_bands[sort_inds[j]]])

        #
        np.fill_diagonal(corrs_mat, 0)
        np.fill_diagonal(ps_mat, 0)


        """ - OLD
        # Theta-Alpha Corr
        [r_th_al, p_th_al] = pearsonr(self.theta_score, self.alpha_score)

        # Theta-Beta Corr
        [r_th_be, p_th_be] = pearsonr(self.theta_score, self.beta_score)

        # Theta-LG Corr
        [r_th_lg, p_th_lg] = pearsonr(self.theta_score, self.lowgamma_score)

        # Alpha-Beta Corr
        [r_al_be, p_al_be] = pearsonr(self.alpha_score, self.beta_score)

        # Alpha-LG Corr
        [r_al_lg, p_al_lg] = pearsonr(self.alpha_score, self.lowgamma_score)

        # Beta-LG Corr
        [r_be_lg, p_be_lg] = pearsonr(self.beta_score, self.lowgamma_score)

        # Save correlation results in a list to return
        corrs = [['Theta-Alpha', r_th_al, p_th_al],
                 ['Theta-Beta' , r_th_be, p_th_be],
                 ['Theta-LG'   , r_th_lg, p_th_lg],
                 ['Alpha-Beta' , r_al_be, p_al_be],
                 ['Alpha-LG'   , r_al_lg, p_al_lg],
                 ['Beta-LG'    , r_be_lg, p_be_lg]]

        # Save corrs out in a matrix
        corrs_mat = np.array([[0, r_th_al, r_th_be, r_th_lg], [r_th_al, 0, r_al_be, r_al_lg],
                              [r_th_be, r_al_be, 0, r_be_lg], [r_th_lg, r_al_lg, r_be_lg, 0]])

        return corrs, corrs_mat
        """


    def save_osc_score(self, file_name):
        """Save out the oscillation score as an npz file.

        Parameters
        ----------
        self : GroupMegData() object.
            Object to store map data across a group of subjects.
        file_name : str
            Name to save the file as.
        """

        # Create full file path and save file as an npz file
        npz_file_name = file_name + '.npz'
        npz_save_name = os.path.join(self.maps_path, 'Oscs', npz_file_name)
        np.savez(npz_save_name, osc_score_theta=self.theta_score, osc_score_alpha=self.alpha_score,
                 osc_score_beta=self.beta_score, osc_score_lowgamma=self.lowgamma_score)


    def set_score_viz(self):
        """Saves a matfile (of oscillation scores) to be loaded for visualization. """

        # Set up paths to save to
        save_name = 'Group_Osc_Score_Viz'
        save_file = os.path.join(self.viz_path, save_name)

        # Save desired outputs into a dictionary
        save_dict = {}
        save_dict['theta_score'] = self.theta_score
        save_dict['alpha_score'] = self.alpha_score
        save_dict['beta_score'] = self.beta_score
        save_dict['lowgamma_score'] = self.lowgamma_score
        save_dict['dat_source'] = self.dat_source
        save_dict['n_subjs'] = self.n_subjs
        save_dict['save_time'] = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

        # Save the dicionary out to a .mat file
        sio.savemat(save_file, save_dict)


    def calc_osc_peak_age(self):
        """Compares age and peak frequency within frequency bands.
        
        Returns
        -------
        corrs : dict
            A dictionary containing correlations results comparing age to oscillations.
        """

        ## NEW

        # Check how many bands there are
        n_bands = len(self.bands)

        # Initialize matrices to store correlation results
        corrs_mat = np.zeros([n_bands])
        ps_mat = np.zeros([n_bands])

        # xx
        sorted_bands, sort_inds = _band_sort(self.bands)

        # xx
        for i in range(n_bands):
            corrs_mat[i], ps_mat[i] = pearsonr(
                self.age, self.peaks[sorted_bands[sort_inds[i]]])

        return corrs_mat, ps_mat, sorted_bands

        """
        # Check correlations
        [r_age_th_peak, p_age_th_peak] = pearsonr(self.age, self.peak_theta)
        [r_age_al_peak, p_age_al_peak] = pearsonr(self.age, self.peak_alpha)
        [r_age_be_peak, p_age_be_peak] = pearsonr(self.age, self.peak_beta)
        [r_age_lg_peak, p_age_lg_peak] = pearsonr(self.age, self.peak_lowgamma)

        # Store correlation results in list to return
        corrs = [['Theta Peak - Age' , r_age_th_peak, p_age_th_peak],
                 ['Alpha Peak - Age' , r_age_al_peak, p_age_al_peak],
                 ['Beta Peak - Age'  , r_age_be_peak, p_age_be_peak],
                 ['LG Peak - Age'    , r_age_lg_peak, p_age_lg_peak]]

        return corrs
        """


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


############################################################################################
########################## OMEGAMAPPIN - OM_MD - PUBLIC FUNCTIONS ##########################
############################################################################################

def print_corrs_mat(rs_mat, ps_mat, labels):
    """Prints out correlations.

    Parameters
    ----------
    rs_mat : 2d array
        xx
    ps_mat : 2d array
        xx
    labels : list(str)
        xx
    """

    # Check how size of the matrix there are
    n = len(labels)

    # ??
    for x in range(n):
        for y in range(n):

            # xx
            if x == y or y < x:
                continue

            print('Corr of ', '{:16}'.format(labels[x]+'-'+labels[y]),
                  ' is ', '{:+1.4f}'.format(rs_mat[x, y]), '    with p-val of ',
                  '{:1.5f}'.format(ps_mat[x, y]))


def print_corrs_vec(rs_vec, ps_vec, labels, desc):
    """   """

    # Check the length of the vector
    n = len(labels)

    #
    for x in range(n):
        print('Corr of ', '{:20}'.format(labels[x]+'-'+desc),' is ',
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


def _get_osc(centers, powers, bws, osc_low, osc_high):
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
    osc_inds = (centers > osc_low) & (centers < osc_high)

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
    dat_source: str
        Which database subject is from.
            Options: {'OMEGA', 'HCP'}

    Returns
    -------
    sex : str
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
    osc_mat : ??
        Matrix of [n_vertex, n_dim, n_subj]

    Returns
    -------
    prob : ??
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
    osc_mat : 3-d array
        XX

    Returns
    -------
    pow_ratio : ??
        xx
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
    avg : str
        What kind of average to take.
            Options: {'mean', 'median'}

    Returns
    -------
    peak : float
        Peak frequency value - the average frequency within a given range.
    """

    #
    osc_inds = (centers > osc_low) & (centers < osc_high)
    osc_cens = centers[osc_inds]

    #
    if avg is 'mean':
        peak = np.mean(osc_cens)
    elif avg is 'median':
        peak = np.median(osc_cens)

    return peak


def _band_sort(osc_bands):
    """

    Parameters
    ----------
    osc_bands : dict
        A dictionary containing the oscillation band definitions.

    Returns
    -------
    ordered_bands : list(str)
        A list of the oscillation band names, in order.
    sort_inds : list(int)
        A list of indices that ...
    """

    # Check how many oscillation bands there are
    n_bands = len(osc_bands)

    # Get low end for each band
    band_names = []
    lower_bounds = np.array([])

    # xx
    for band in osc_bands:
        band_names.append(band)
        lower_bounds = np.append(lower_bounds, osc_bands[band][0])

    #
    sort_inds = np.argsort(lower_bounds)

    #
    ordered_bands = []
    ordered_bands[:] = [band_names[i] for i in sort_inds]

    return ordered_bands, sort_inds


def _load_foof_pickle(path):
    """Loads FOOF data from a pickle file.

    Parameters
    ----------
    path : str
        xx
    """

    results = pickle.load(open(path, 'rb'))
    return results


def _load_foof_csv(path):
    """
    NOTE: not yet implemented
    """

    pass
