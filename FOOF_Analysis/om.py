from __future__ import print_function, division
from scipy.stats.stats import pearsonr
import numpy as np
import os
import csv
import scipy.io as sio
import pickle
import matplotlib.pyplot as plt


########################################################################################
###############################  OMEGAMAPPING - CLASSES  ###############################
########################################################################################

class Osc:
	""" Class to store oscillations parameters. 
	"""
	def __init__(self):

		# Theta
		self.theta_low = 3
		self.theta_high = 8
		
		# Alpha
		self.alpha_low = 8
		self.alpha_high = 13
		
		# Beta
		self.beta_low = 13
		self.beta_high = 30

		# Low Gamma
		self.lowgamma_low = 30
		self.lowgamma_high = 40


class MegData():
	""" Class for a single subject of FOOF results for MEG Source PSDs. 
	"""

	def __init__(self):

		# Initialize subject number
		self.subnum = int()
		self.nPSDs = int()

		# Initialize data vectors
		self.chis = np.array([])
		self.centers = np.array([])
		self.powers = np.array([])
		self.bws = np.array([])

		# Initialize vectors for all-osc data
		self.centers_all = np.array([])
		self.powers_all = np.array([])
		self.bws_all = np.array([])
		self.nOscs = np.array([])

		# Initialize matrices for osc-band data
		self.gr_thetas = np.array([])
		self.gr_alphas = np.array([])
		self.gr_betas = np.array([])
		self.gr_lowgammas = np.array([])

		# Set plot title
		self.title = ''
		self.vis_opac = 0.1

		# Set boolean for what has been run
		self.has_data = False
		self.all_osc = False
		self.bands_vertex = False

		# Initialize demographic variables
		self.sex = []
		self.age = np.array([])

		# Initialize peak frequency variables
		self.peak_theta = np.array([])
		self.peak_alpha = np.array([])
		self.peak_beta = np.array([])
		self.peak_lowgamma = np.array([])


	def import_foof(self, subnum, get_demo=True):
		""" Import FOOF results from pickle file. 
		"""

		# Check if object already has data
		if(self.has_data):
			print("Subject object already contains add. Can't add")
			return

		# Set subject number for current data object
		self.subnum = subnum
		self.title = 'S-' + str(self.subnum)

		#
		foof_data_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/3-FOOF/'
		files = os.listdir(foof_data_path)
		files = clean_file_list(files, 'Foof_Vertex')

		#
		cur_subj_file = get_cur_subj(subnum, files)
		cur_subj_path = os.path.join(foof_data_path, cur_subj_file)

		# Load pickle file and check the number of PSDs
		results = pickle.load(open(cur_subj_path, 'rb'))
		self.nPSDs = len(results)

		# Initialize numpy arrays to pull out different result params
		self.chis = np.zeros([self.nPSDs, 1])
		self.centers = np.zeros([self.nPSDs, 8])
		self.powers = np.zeros([self.nPSDs, 8])
		self.bws = np.zeros([self.nPSDs, 8])

		# Loop through FOOF results, pulling out individual findings
		for i in range(0, self.nPSDs):
		    self.chis[i, 0] = results[i][0]
		    self.centers[i, 0:len(results[i][1])] = results[i][1]
		    self.powers[i, 0:len(results[i][2])] = results [i][2]
		    self.bws[i, 0:len(results[i][3])] = results[i][3]

		# Check how many oscillations per vertex
		self.osc_count = np.zeros([self.nPSDs, 1])
		for i in range(0, self.nPSDs):
			self.osc_count[i, 0] = len(np.nonzero(self.centers[i, :])[0])

		# Get demographic data
		if(get_demo):
			self.sex, self.age = _get_demo_csv(self.subnum)

		# Update boolean to say current subject has data attached
		self.has_data = True


	def osc_bands_vertex(self, osc):
		"""
		"""

		## Re-Initialize matrices to right size to save results
		self.thetas = np.zeros([self.nPSDs, 4])
		self.alphas = np.zeros([self.nPSDs, 4])
		self.betas = np.zeros([self.nPSDs, 4])
		self.lowgammas = np.zeros([self.nPSDs, 4])

		#
		for i in range(0, self.nPSDs):

			#
			centers_temp = self.centers[i, :]
			powers_temp = self.powers[i, :]
			bws_temp = self.bws[i, :]

			self.thetas[i, :] = _get_osc(centers_temp, powers_temp, bws_temp, 
				osc.theta_low, osc.theta_high)

			self.alphas[i, :] = _get_osc(centers_temp, powers_temp, bws_temp, 
				osc.alpha_low, osc.alpha_high)

			self.betas[i, :] = _get_osc(centers_temp, powers_temp, bws_temp, 
				osc.beta_low, osc.beta_high)

			self.lowgammas[i, :] = _get_osc(centers_temp, powers_temp, bws_temp, 
				osc.lowgamma_low, osc.lowgamma_high)

		self.bands_vertex = True


	def save_viz(self):
		""" Saves a matfile of frequency information to be loaded with Brainstorm for visualization. 
		"""

		#
		save_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/4-Viz/'
		save_name = str(self.subnum) + '_Foof_Viz'
		save_file = os.path.join(save_path, save_name)

		#
		save_dict = {}
		save_dict['slopes'] = self.chis
		save_dict['thetas'] = self.thetas
		save_dict['alphas'] = self.alphas
		save_dict['betas'] = self.betas
		save_dict['lowgammas'] = self.lowgammas

		#
		sio.savemat(save_file, save_dict)


	def all_oscs(self):
		""" Flatten osc data to vectors. 
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

		# Get the number of oscillations
		self.nOscs = len(self.centers_all)

		# Update format
		self.all_osc = True


	def peak_freq(self, osc):

		#
		self.peak_theta = _osc_peak(self.centers_all, osc.theta_low, osc.theta_high)
		self.peak_alpha = _osc_peak(self.centers_all, osc.alpha_low, osc.alpha_high)
		self.peak_beta = _osc_peak(self.centers_all, osc.beta_low, osc.beta_high)
		self.peak_lowgamma = _osc_peak(self.centers_all, osc.lowgamma_low, osc.lowgamma_high)


	def plot_all_oscs(self):
		""" Plots histogram distributions of oscillation centers, powers and bws. 
		"""

		# Check if all_oscs computed
		if(not self.all_osc):
			print("All oscillations not computed. Can't plot.")
			return

		# Plot Settings
		nBins = 160 	# Number of bins for histograms
		st_fs = 20 		# Super Title Font Size
		sp_fs = 18		# Subplit Title Font Size
		ax_fs = 16		# Axis Label Font Size

		# Set up subplots
		f, ax = plt.subplots(3, 1, figsize=(15, 15))
		# Set plot super-title
		plt.suptitle('Distributions of Oscillatory Parameters - ' + self.title, fontsize=st_fs, fontweight='bold')

		# Subplot 1 - Center Frequency
		ax[0].hist(self.centers_all, nBins)
		ax[0].set_title('Center Frequency', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[0].set_xlabel('Frequency', {'fontsize': ax_fs})
		ax[0].set_ylabel('Count', {'fontsize': ax_fs})

		# Subplot 2 - Power
		ax[1].hist(np.log10(self.powers_all), nBins)
		ax[1].set_title('Oscillatory Power', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[1].set_xlabel('Log Power', {'fontsize': ax_fs})
		ax[1].set_ylabel('Count', {'fontsize': ax_fs})

		# Subplot 3 - Bandwidth
		ax[2].hist(self.bws_all, nBins)
		ax[2].set_title('Band Width', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[2].set_xlabel('Bandwidth (Hz)', {'fontsize': ax_fs})
		ax[2].set_ylabel('Count', {'fontsize': ax_fs})

		# Adjust subplot spacing
		plt.subplots_adjust(hspace=0.4)


	def plot_hist_count(self):
		""" Plots a histogram of the osc_count vector. 
		"""

		# Plot Settings
		nBins = 25
		t_fs = 18
		ax_fs = 16

		# Create histogram
		plt.hist(self.osc_count, nBins, range=[0, 8])
		plt.title('# Oscillations per Vertex', {'fontsize': t_fs, 'fontweight': 'bold'})
		plt.xlabel('# Oscillations', {'fontsize': ax_fs, 'fontweight': 'bold'})
		plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})


	def plot_slopes(self):
		""" Plots a histogram of the chi values for all vertices. 
		"""

		# Plot Settings
		nBins = 100
		t_fs = 20
		ax_fs = 16

		# Create histogram
		plt.hist(self.chis, nBins)
		plt.title('Slopes - ' + self.title, {'fontsize': t_fs, 'fontweight': 'bold'})
		plt.xlabel('Chi Parameter', {'fontsize': ax_fs, 'fontweight': 'bold'})
		plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})


	def plot_comparison(self):

		# Check if all-osc computed
		if(not self.all_osc):
			print("All oscillations not computed. Can't run comparisons.")
			return

		# Plot Settings
		st_fs = 20 		# Super Title Font Size
		sp_fs = 18		# Subplit Title Font Size
		ax_fs = 16		# Axis Label Font Size

		# Set up subplots
		f, ax = plt.subplots(3, 1, figsize=(15, 15))
		# Set plot super-title
		plt.suptitle('Oscillation Parameter Comparisons - ' + self.title, fontsize=st_fs, fontweight='bold')

		## Centers vs. Bandwidth
		# Check Correlation
		corr_cen_bw, pcorr_cen_bw = pearsonr(self.centers_all, np.log10(self.bws_all))

		# Plot
		ax[0].plot(self.centers_all, np.log10(self.bws_all), '.', alpha=self.vis_opac)
		ax[0].set_title('Center vs. Bandwidth', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[0].set_xlabel('Centers', {'fontsize': ax_fs})
		ax[0].set_ylabel('BW', {'fontsize': ax_fs})


		## Centers vs. Power
		# Check Correlation
		corr_cen_pow, pcorr_cen_pow = pearsonr(self.centers_all, np.log10(self.powers_all))

		# Plot
		ax[1].plot(self.centers_all, np.log10(self.powers_all), '.', alpha=self.vis_opac)
		ax[1].set_title('Center vs. Power', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[1].set_xlabel('Centers', {'fontsize': ax_fs})
		ax[1].set_ylabel('Log Power', {'fontsize': ax_fs})


		## Bandwidth vs. Power
		# Check Correlation
		corr_bw_pow, pcorr_bw_pow = pearsonr(np.log10(self.bws_all), np.log10(self.powers_all))

		# Plot
		ax[2].plot(np.log10(self.bws_all), np.log10(self.powers_all), '.', alpha=self.vis_opac)		
		ax[2].set_title('BW vs. Power', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[2].set_xlabel('Bandwidth (Hz)', {'fontsize': ax_fs})
		ax[2].set_ylabel('Log Power', {'fontsize': ax_fs})

		# Adjust subplot spacing
		plt.subplots_adjust(hspace=0.4)

		corrs = [['Center - B.W.'   , corr_cen_bw , pcorr_cen_bw ],
				 ['Center - Power'  , corr_cen_pow, pcorr_cen_pow],
				 ['B.W. - Power '   , corr_bw_pow , pcorr_bw_pow ],
				 ]

		return corrs


	def save_pickle(self):
		""" Save current meg data object as a pickled object. 
		"""
		pass


class GroupMegData(MegData):
	""" A class to store OMEGA data from multiple subjects.
	Holds all oscillations, regardless of spatial location. 

	Note: Class derived from MegData
	"""


	def __init__(self):

		MegData.__init__(self)

		# Initialize groups subject variables
		self.nsubjs = int()
		self.subjs = []

		# 
		self.nOscs_tot = int()

		# Set title for plots
		self.title = 'Group' 
		self.vis_opac = 0.005

		# Set booleans for what has been run
		self.osc_prob_done = False
		self.osc_score_done = False


	def add_subject(self, new_subj, add_all_oscs=False, add_vertex_bands=False, add_vertex_oscs=False):
		"""
		Inputs:
			new_subj	- MEG subject (instance of MegData)

		"""

		# Check if subject has data
		if(not new_subj.has_data):
			print("Empty meg data object. Cannot add data.")

		# Add All-Osc Data
		if(add_all_oscs):
			
			# Add oscillation parameters to current data 
			self.centers_all = np.append(self.centers_all, new_subj.centers_all)
			self.bws_all = np.append(self.bws_all, new_subj.bws_all)
			self.powers_all = np.append(self.powers_all, new_subj.powers_all)
			self.chis = np.append(self.chis, new_subj.chis)

			# Update count of total number of oscillations
			self.nOscs = np.append(self.nOscs, new_subj.nOscs)
			self.nOscs_tot = len(self.centers_all)

		# 
		if(add_vertex_bands):
			if(self.nsubjs == 0):
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
		if(add_vertex_oscs):
			if(self.nsubjs == 0):
				self.centers = new_subj.centers
				self.powers = new_subj.powers
				self.bws = new_subj.bws
			else:
				self.centers = np.dstack([self.centers, new_subj.centers])
				self.powers = np.dstack([self.powers, new_subj.powers])
				self.bws = np.dstack([self.bws, new_subj.bws])

		# Update subj count and subject number list
		self.nsubjs += 1
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
		"""
		"""

		#
		self.theta_prob = _osc_prob(self.gr_thetas)
		self.alpha_prob = _osc_prob(self.gr_alphas)
		self.beta_prob = _osc_prob(self.gr_betas)
		self.lowgamma_prob = _osc_prob(self.gr_lowgammas)

		#
		self.osc_prob_done = True


	def set_prob_vis(self):
		""" Saves a matfile (of osc-probs) to be loaded with Brainstorm for visualization. 
		"""

		#
		save_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/4-Viz/'
		save_name = 'Group_Osc_Prob_Viz'
		save_file = os.path.join(save_path, save_name)

		#
		save_dict = {}
		save_dict['theta_prob'] = self.theta_prob
		save_dict['alpha_prob'] = self.alpha_prob
		save_dict['beta_prob'] = self.beta_prob
		save_dict['lowgamma_prob'] = self.lowgamma_prob

		#
		sio.savemat(save_file, save_dict)


	def osc_prob_corrs(self):
		"""
		"""

		# Check if oscillation probabilities have been calculated. 
		if(not self.osc_prob_done):
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

		corrs = [['Theta-Alpha', r_th_al, p_th_al],
				 ['Theta-Beta' , r_th_be, p_th_be],
				 ['Theta-LG'   , r_th_lg, p_th_lg],
				 ['Alpha-Beta' , r_al_be, p_al_be],
				 ['Alpha-LG'   , r_al_lg, p_al_lg],
				 ['Beta-LG'    , r_be_lg, p_be_lg],
				 ]

		return corrs


	def osc_score_corrs(self):
		"""
		"""

		# Check if oscillation probabilities have been calculated. 
		if(not self.osc_score_done):
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

		corrs = [['Theta-Alpha', r_th_al, p_th_al],
				 ['Theta-Beta' , r_th_be, p_th_be],
				 ['Theta-LG'   , r_th_lg, p_th_lg],
				 ['Alpha-Beta' , r_al_be, p_al_be],
				 ['Alpha-LG'   , r_al_lg, p_al_lg],
				 ['Beta-LG'    , r_be_lg, p_be_lg],
				 ]

		return corrs


	def osc_score(self):

		#
		if(not self.osc_prob_done):
			print("Oscillation probability not computed. Can't continue.")

		#	
		self.theta_pow_ratio = _osc_pow_ratio(self.gr_thetas)
		self.alpha_pow_ratio = _osc_pow_ratio(self.gr_alphas)
		self.beta_pow_ratio = _osc_pow_ratio(self.gr_betas)
		self.lowgamma_pow_ratio = _osc_pow_ratio(self.gr_lowgammas)

		#
		self.theta_score = self.theta_pow_ratio * self.theta_prob
		self.alpha_score = self.alpha_pow_ratio * self.alpha_prob
		self.beta_score = self.beta_pow_ratio * self.beta_prob
		self.lowgamma_score = self.lowgamma_pow_ratio * self.lowgamma_prob

		#
		self.osc_score_done = True


	def set_score_vis(self):
		""" Saves a matfile (of osc-score) to be loaded with Brainstorm for visualization. 
		"""

		#
		save_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/4-Viz/'
		save_name = 'Group_Osc_Score_Viz'
		save_file = os.path.join(save_path, save_name)

		#
		save_dict = {}
		save_dict['theta_score'] = self.theta_score
		save_dict['alpha_score'] = self.alpha_score
		save_dict['beta_score'] = self.beta_score
		save_dict['lowgamma_score'] = self.lowgamma_score

		#
		sio.savemat(save_file, save_dict)


	def osc_age_comparison_plot(self):

		# Plot settings
		st_fs = 20
		sp_fs = 18

		# Set up subplots
		f, ax = plt.subplots(2, 2, figsize=(10, 10))
		# Set plot super-title
		plt.suptitle('Peak Frequency / Age Comparisons', fontsize=st_fs, fontweight='bold')

		# Theta
		ax[0,0].plot(self.age, self.peak_theta, '.')
		ax[0,0].set_title('Theta', {'fontsize': sp_fs, 'fontweight': 'bold'})

		# Alpha
		ax[0,1].plot(self.age, self.peak_alpha, '.')
		ax[0,1].set_title('Alpha', {'fontsize': sp_fs, 'fontweight': 'bold'})

		# Beta
		ax[1,0].plot(self.age, self.peak_beta, '.')
		ax[1,0].set_title('Beta', {'fontsize': sp_fs, 'fontweight': 'bold'})

		# Gamma
		ax[1,1].plot(self.age, self.peak_lowgamma, '.')
		ax[1,1].set_title('Low Gamma', {'fontsize': sp_fs, 'fontweight': 'bold'})

		# Check correlations
		[r_age_th_peak, p_age_th_peak] = pearsonr(self.age, self.peak_theta)
		[r_age_al_peak, p_age_al_peak] = pearsonr(self.age, self.peak_alpha)
		[r_age_be_peak, p_age_be_peak] = pearsonr(self.age, self.peak_beta)
		[r_age_lg_peak, p_age_lg_peak] = pearsonr(self.age, self.peak_lowgamma)

		# Store correlation results in list to return
		corrs = [['Theta Peak - Age' , r_age_th_peak, p_age_th_peak],
				 ['Alpha Peak - Age' , r_age_al_peak, p_age_al_peak],
				 ['Beta Peak - Age'  , r_age_be_peak, p_age_be_peak],
				 ['LG Peak - Age'    , r_age_lg_peak, p_age_lg_peak],
				 ]

		return corrs


	def freq_corr(self, f_win):
		"""

		To test more explicitly, you can take all neighboring pairs of frequencies using, say, 3 Hz segments, and look at correlation coefficient.

		So 2-5 Hz v 5-8 Hz, then 3-6 v 6-9, then 4-7 v 8-11, and so on, all the way up to 33-36 v 37-40.
		"""

		#
		[nVertex, nSlots, nSubj] = np.shape(self.centers)

		nFs = len(range(3, 40-f_win))

		prob_mat = np.zeros([nVertex, nFs])

		#
		for v in range(0, nVertex):

			#
			for s in range(0, nSubj):

				cens_temp = self.centers[v, :, s]

				# 
				i  = 0
				for f in range(3, 40-f_win):

					#
					cens_fwin = _get_all_osc(cens_temp, f, f + f_win)

					#
					if(len(cens_fwin) != 0):
						prob_mat[v, i] += 1

					i += 1

		#
		prob_mat = prob_mat/nSubj

		# 
		corr_vec = np.zeros([nFs-1])
		p_vec = np.zeros([nFs-1])

		# Compute corr between f and f+1 start windows
		#for f in range(0, nFs-1):
			#corr_vec[f], p_vec[f] = pearsonr(prob_mat[:, f], prob_mat[:, f+1])

		# Compute corr between f and f+f_win start windows
		for f_ind in range(0, nFs-f_win):
			corr_vec[f_ind], p_vec[f_ind] = pearsonr(prob_mat[:, f_ind], prob_mat[:, f_ind+f_win])

		return corr_vec, p_vec


########################################################################################
############################### OMEGAMAPPING - FUNCTIONS ###############################
########################################################################################


def clean_file_list(files_in, string):
    '''Takes a list of file names (strings), returns only those with 'string' in them.
    '''
    
    files_out = []

    #
    for i in range(0, len(files_in)):
        if(string in files_in[i]):
            files_out.append(files_in[i])
            
    return files_out


def get_sub_nums(files_in):
	""" Takes a list of files. Returns a list of subject numbers. 
	"""

	sub_nums = []

	#
	for f in files_in:
		str_split = f.split('_', 1)
		sub_nums.append(int(str_split[0]))

	return sub_nums


def get_cur_subj(cur_subj, files):
    '''Takes an int, and a list of file names (strings), returns the file name with the number in it.
    '''
    
    #
    cur_subj_str = str(cur_subj)
    
    #
    for i in range(0, len(files)):
        if(cur_subj_str in files[i]):
            return files[i]


def run_par_foof():
	""" NOT YET IMPLEMENTED.
	"""

	pass


########################################################################################
###############################  OM - PRIVATE FUNCTIONS  ###############################
########################################################################################


def _get_osc(centers, powers, bws, osc_low, osc_high):
	""" Searches for an oscillations of specified frequency band. 
	Returns a single oscillation in that band. 
	Helper function for osc_per_vertex in MegData.
	"""

	#
	osc_inds = (centers > osc_low) & (centers < osc_high)

	#
	osc_cens = centers[osc_inds]
	osc_pows = powers[osc_inds]
	osc_bws = bws[osc_inds]

	# Get number of oscillations in the frequency band. 
	n = len(osc_cens)

	# OLD: Get single oscillation. This just returned the first (lowest freq) osc.
	#cen, n = _get_single_osc(osc_cens)
	#power, x = _get_single_osc(osc_pows)
	#bw, x = _get_single_osc(osc_bws)

	# Get highest power oscillation in band
	cen, power, bw = _get_single_osc_power(osc_cens, osc_pows, osc_bws)

	return np.array([cen, power, bw, n])


def _get_all_osc(centers, osc_low, osc_high):
	"""
	"""

	osc_inds = (centers > osc_low) & (centers < osc_high)
	osc_cens = centers[osc_inds]

	return osc_cens


def _get_demo_csv(sub_num):
	"""
	"""

	#
	csv_data_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/'
	csv_file_name = '00-Collin_Subjects.csv'
	csv_file = os.path.join(csv_data_path, csv_file_name)

	#
	with open(csv_file, 'rb') as f:
		reader = csv.reader(f, delimiter=',')
		for row in reader:
			if(row[1] == str(sub_num)):
				sex = row[4]
				age = int(row[7])
				break

	return sex, age


def _get_single_osc(osc_in):
	''' OLD - UNUSED
	Takes a vector, returns the first element, and the length.
	This selects the lowest frequency oscillation. 
	Helper function for osc_per_vertex in MegData.
	'''

	#
	if(len(osc_in) == 0):
		return 0, 0
	elif(len(osc_in) == 1):
		return osc_in, 1
	else:
		return osc_in[0], len(osc_in)


def _get_single_osc_power(osc_cens, osc_pows, osc_bws):
	"""
	"""

	# 
	if(len(osc_cens) == 0):
		return 0., 0., 0.
	elif(len(osc_cens) == 1):
		return osc_cens, osc_pows, osc_bws
	else:
		high_ind = np.argmax(osc_pows)
		return osc_cens[high_ind], osc_pows[high_ind], osc_bws[high_ind]


def _osc_prob(osc_mat):
	""" Takes a 3D matrix of oscillations across subjects.
	osc_mat: [nVertex, nDim, nSubj]
	Returns a vector of probability of oscillation at each vertex.
	"""

	# 
	[nVertex, nDim, nSubj] = np.shape(osc_mat)

	#
	prob = np.zeros([nVertex])

	#
	for i in range(0, nVertex):
		prob[i] = (np.count_nonzero(osc_mat[i, 0, :]) / nSubj)

	return prob


def _osc_pow_ratio(osc_mat):
	"""
	"""

	#
	[nVertex, nDim, nSubj] = np.shape(osc_mat)

	#
	avg_powers = np.zeros(nVertex)
	for v in range(0, nVertex):
		temp_pows = osc_mat[v, 1, :]
		temp_pows = temp_pows[np.nonzero(temp_pows)]
		if(len(temp_pows) == 0):
			avg_powers[v] = 0
		else:
			avg_powers[v] = np.mean(temp_pows)
	max_all = max(avg_powers)

	#
	pow_ratio = np.zeros(nVertex)
	for v in range(0, nVertex):
		pow_ratio[v] = np.mean(osc_mat[v, 1, :]) / max_all

	return pow_ratio


def _osc_peak(centers, osc_low, osc_high):
	"""
	"""

	#
	osc_inds = (centers > osc_low) & (centers < osc_high)
	osc_cens = centers[osc_inds]
	peak = np.mean(osc_cens)

	return peak
