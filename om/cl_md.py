from __future__ import print_function, division
import os
import csv
import pickle
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
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

        # Initialize arrays for oscillation bands
        self.thetas = np.array([])
        self.alphas = np.array([])
        self.betas = np.array([])
        self.lowgammas = np.array([])

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

        # Initialize peak frequency variables
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


    def osc_bands_vertex(self, osc):
        """Groups oscillations at each vertex in distinct frequency bands.
        Stores band specific oscillations in (self.){thetas, alphas, betas, lowgammas}.

        Parameters
        ----------
        self : MegData object
            MegData object.
        osc : Osc object
            An object containing frequency bands to use.
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


    def save_viz(self):
        """Saves a matfile of freq info to be loaded with Brainstorm for visualization. """

        # Set up paths to save to
        #save_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/4-Viz/'
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


    def all_oscs(self):
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
            print('Subj:', str(self.subnum), 'Found', str(n_nans), ' NaNs. Removing.')

            # Remove osc's with nan values
            non_nans = np.logical_not(nans)
            self.centers_all = self.centers_all[non_nans]
            self.powers_all = self.powers_all[non_nans]
            self.bws_all = self.bws_all[non_nans]

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

        # Get peak frequency within each frequency band
        self.peak_theta = _osc_peak(self.centers_all, osc.theta_low, osc.theta_high, avg)
        self.peak_alpha = _osc_peak(self.centers_all, osc.alpha_low, osc.alpha_high, avg)
        self.peak_beta = _osc_peak(self.centers_all, osc.beta_low, osc.beta_high, avg)
        self.peak_lowgamma = _osc_peak(self.centers_all, osc.lowgamma_low, osc.lowgamma_high, avg)


    def calc_osc_param_corrs(self):
        """Calculates correlations between oscillatory parameters."""

        ## Calculate correlations

        # Centers vs. Bandwidth
        corr_cen_bw, pcorr_cen_bw = pearsonr(self.centers_all, np.log10(self.bws_all))

        # Centers vs. Power
        corr_cen_pow, pcorr_cen_pow = pearsonr(self.centers_all, np.log10(self.powers_all))

        # Bandwidth vs. Power
        corr_bw_pow, pcorr_bw_pow = pearsonr(np.log10(self.bws_all), np.log10(self.powers_all))

        # Save correlation results to list to return
        corrs = [['Center - B.W.'   , corr_cen_bw , pcorr_cen_bw ],
                 ['Center - Power'  , corr_cen_pow, pcorr_cen_pow],
                 ['B.W. - Power '   , corr_bw_pow , pcorr_bw_pow ]]

        return corrs


    def save_pickle(self):
        """Save current meg data object as a pickled object.

        NOTE: NOT YET IMPLEMENTED.
        """

        pass


class GroupMegData(MegData):
    """A class to store OMEGA data from multiple subjects.

    Holds all oscillations, regardless of spatial location.
    Note: Class derived from MegData
    """

    def __init__(self, db):

        #
        MegData.__init__(self, db)

        # Initialize groups subject variables
        self.n_subjs = int()
        self.subjs = []

        #
        self.n_oscs_tot = int()

        # Set title for plots
        self.title = 'Group'

        # Initialize matrices for osc-band data
        self.gr_thetas = np.array([])
        self.gr_alphas = np.array([])
        self.gr_betas = np.array([])
        self.gr_lowgammas = np.array([])

        # Initialize to store oscillation probabilities
        self.theta_prob = np.array([])
        self.alpha_prob = np.array([])
        self.beta_prob = np.array([])
        self.lowgamma_prob = np.array([])

        # Initialize to store oscillation power ratios
        self.theta_pow_ratio = np.array([])
        self.alpha_pow_ratio = np.array([])
        self.beta_pow_ratio = np.array([])
        self.lowgamma_pow_ratio = np.array([])

        # Initialize to store oscillation scores
        self.theta_score = np.array([])
        self.alpha_score = np.array([])
        self.beta_score = np.array([])
        self.lowgamma_score = np.array([])

        # Initialize vars to store slope values
        self.gr_chis = np.array([])
        self.slopes_avg = np.array([])

        # Set booleans for what has been run
        self.osc_prob_done = False
        self.osc_score_done = False


    def add_subject(self, new_subj, add_all_oscs=False, add_vertex_bands=False, 
                    add_vertex_oscs=False, add_vertex_slopes=False):
        """Adds a new subject to the GroupMegData object.

        Parameters
        ----------
        self : GroupMegData() object. 
            Object to store map data across a group of subjects.
        new_subj : MegData() Object
            MEG subject (instance of MegData)
        add_all_oscs : boolean, optional
            XX
        add_vertex_bands : boolean, optional
            XX
        add_vertex_oscs : boolean, optional
            XX
        add_vertex_slopes : boolean, optional
            xx
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

            # Update count of total number of oscillations
            self.n_oscs = np.append(self.n_oscs, new_subj.n_oscs)
            self.n_oscs_tot = len(self.centers_all)

        # Add band-specific data
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

        #
        if add_vertex_slopes:

            if self.n_subjs == 0:
                self.gr_chis = new_subj.slopes
            else:
                self.gr_chis = np.hstack([self.gr_chis, new_subj.slopes])

        # Update subj count and subject number list
        self.n_subjs += 1
        self.subjs = np.append(self.subjs, new_subj.subnum)

        #
        self.all_osc = add_all_oscs
        self.bands_vertex = add_vertex_bands

        # Add demographic data
        self.sex.append(new_subj.sex)
        self.age = np.append(self.age, new_subj.age)

        # Add osc-peak data
        self.peak_theta = np.append(self.peak_theta, new_subj.peak_theta)
        self.peak_alpha = np.append(self.peak_alpha, new_subj.peak_alpha)
        self.peak_beta = np.append(self.peak_beta, new_subj.peak_beta)
        self.peak_lowgamma = np.append(self.peak_lowgamma, new_subj.peak_lowgamma)


    def osc_prob(self):
        """Calculates the probability (per vertex / across subjects) of an osc in a specific band. """

        # For each frequency band, compute the probability of oscillation in that band.
        self.theta_prob = _osc_prob(self.gr_thetas)
        self.alpha_prob = _osc_prob(self.gr_alphas)
        self.beta_prob = _osc_prob(self.gr_betas)
        self.lowgamma_prob = _osc_prob(self.gr_lowgammas)

        # Update boolean that oscillation probability has been computed
        self.osc_prob_done = True


    def group_slope(self, save_out=False, file_name=None, set_viz=False):
        """Calculates the average slope value for each vertex, across subjects.

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
        self.slopes_avg = np.median(self.gr_chis, axis=1)

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
            save_dict['n_subjs'] = self.n_subjs

            # Save the dicionary out to a .mat file
            sio.savemat(save_file, save_dict)


    def set_prob_vis(self):
        """Saves out a matfile (of osc-probs) to be loaded with Brainstorm for visualization. """

        # Set up paths to save to
        save_name = 'Group_Osc_Prob_Viz'
        save_file = os.path.join(self.viz_path, save_name)

        # Save desired outputs into a dictionary
        save_dict = {}
        save_dict['theta_prob'] = self.theta_prob
        save_dict['alpha_prob'] = self.alpha_prob
        save_dict['beta_prob'] = self.beta_prob
        save_dict['lowgamma_prob'] = self.lowgamma_prob

        # Save the dicionary out to a .mat file
        sio.savemat(save_file, save_dict)


    def osc_prob_corrs(self):
        """Calculates the correlations between oscillation probabilities.

        Returns
        -------
        corrs : dict
            Contains the correlations between all oscillation probabilities.
        """

        # Check if oscillation probabilities have been calculated.
        if not self.osc_prob_done:
            print("Osc Probability not calculated - can't proceed.")
            return

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

        return corrs


    def osc_score(self):
        """Calculate the oscillation score for each frequency band.

        The oscillation score is .... XXXXX

        """

        # Check if osc-prob is calculated. Can't proceed if it isnt.
        if not self.osc_prob_done:
            print("Oscillation probability not computed. Can't continue.")

        # Compute power ratio for each oscillatory band
        self.theta_pow_ratio = _osc_pow_ratio(self.gr_thetas)
        self.alpha_pow_ratio = _osc_pow_ratio(self.gr_alphas)
        self.beta_pow_ratio = _osc_pow_ratio(self.gr_betas)
        self.lowgamma_pow_ratio = _osc_pow_ratio(self.gr_lowgammas)

        # Compute oscillation score for each oscillatory band
        self.theta_score = self.theta_pow_ratio * self.theta_prob
        self.alpha_score = self.alpha_pow_ratio * self.alpha_prob
        self.beta_score = self.beta_pow_ratio * self.beta_prob
        self.lowgamma_score = self.lowgamma_pow_ratio * self.lowgamma_prob

        # Set boolean that osc-score has been computed
        self.osc_score_done = True


    def osc_score_corrs(self):
        """Calculates the correlations between oscillations scores.

        Returns
        -------
        corrs : dict
            Contains the results of the correlations across oscillation scores.
        """

        # Check if oscillation probabilities have been calculated.
        if not self.osc_score_done:
            print("Osc Score not calculated - can't proceed.")
            return

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

        return corrs


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


    def set_score_vis(self):
        """Saves a matfile (of osc-scores) to be loaded with Brainstorm for visualization. """

        # Set up paths to save to
        save_name = 'Group_Osc_Score_Viz'
        save_file = os.path.join(self.viz_path, save_name)

        # Save desired outputs into a dictionary
        save_dict = {}
        save_dict['theta_score'] = self.theta_score
        save_dict['alpha_score'] = self.alpha_score
        save_dict['beta_score'] = self.beta_score
        save_dict['lowgamma_score'] = self.lowgamma_score

        # Save the dicionary out to a .mat file
        sio.savemat(save_file, save_dict)


    def calc_osc_peak_age(self):
        """Compares age and peak frequency within frequency bands.
        
        Returns
        -------
        corrs : dict
            A dictionary containing correlations results comparing age to oscillations.

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

    Power ratio is a score between [0, 1] power relative to
        max power in that frequency band.

    Parameters
    ----------
    osc_mat : ??
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
