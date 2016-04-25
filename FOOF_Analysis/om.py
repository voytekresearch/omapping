from __future__ import print_function, division
from scipy.stats.stats import pearsonr
import numpy as np
import os
import scipy.io as sio
import pickle
import matplotlib.pyplot as plt


########################################################################
#######################  OMEGAMAPPING - CLASSES  #######################
########################################################################

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

		# Initialize matrices for osc-band data
		self.gr_thetas = np.array([])
		self.gr_alphas = np.array([])
		self.gr_betas = np.array([])
		self.gr_lowgammas = np.array([])

		# Set plot title
		self.title = ''

		# Set boolean for what has been run
		self.has_data = False
		self.all_osc = False
		self.bands_vertex = False

	def import_foof(self, subnum):
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

		# 
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
		""" Saves a matfile in matfile format in Brainstorm for visualization. 
		"""
		save_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/4-Viz/'
		save_name = str(self.subnum) + '_Foof_Viz'
		save_file = os.path.join(save_path, save_name)

		save_dict = {}
		save_dict['slopes'] = self.chis
		save_dict['thetas'] = self.thetas
		save_dict['alphas'] = self.alphas
		save_dict['betas'] = self.betas
		save_dict['lowgammas'] = self.lowgammas

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
		self.nOscs = len(self.centers)

		# Update format
		self.all_osc = True

	def plot_all_oscs(self):
		""" Plots histogram distributions of oscillation centers, powers and bws. 
		"""

		# Check if all_oscs computed
		if(not self.all_osc):
			print("All oscillations not computed. Can't plot.")
			return

		# Plot Settings
		nBins = 100 	# Number of bins for histograms
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
		vis_alp = 0.1	# Opacity for plotting

		# Set up subplots
		f, ax = plt.subplots(3, 1, figsize=(15, 15))
		# Set plot super-title
		plt.suptitle('Oscillation Parameter Comparisons - ' + self.title, fontsize=st_fs, fontweight='bold')

		## Centers vs. Bandwidth
		# Check Correlation
		corr_cen_bw, pcorr_cen_bw = pearsonr(self.centers_all, np.log10(self.bws_all))

		# Plot
		ax[0].plot(self.centers_all, np.log10(self.bws_all), '.', alpha=vis_alp)
		ax[0].set_title('Center vs. Bandwidth', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[0].set_xlabel('Centers', {'fontsize': ax_fs})
		ax[0].set_ylabel('BW', {'fontsize': ax_fs})


		## Centers vs. Power
		# Check Correlation
		corr_cen_pow, pcorr_cen_pow = pearsonr(self.centers_all, np.log10(self.powers_all))

		# Plot
		ax[1].plot(self.centers_all, np.log10(self.powers_all), '.', alpha=vis_alp)
		ax[1].set_title('Center vs. Power', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[1].set_xlabel('Centers', {'fontsize': ax_fs})
		ax[1].set_ylabel('Log Power', {'fontsize': ax_fs})


		## Bandwidth vs. Power
		# Check Correlation
		corr_bw_pow, pcorr_bw_pow = pearsonr(np.log10(self.bws_all), np.log10(self.powers_all))

		# Plot
		ax[2].plot(np.log10(self.bws_all), np.log10(self.powers_all), '.', alpha=vis_alp)		
		ax[2].set_title('BW vs. Power', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[2].set_xlabel('Bandwidth (Hz)', {'fontsize': ax_fs})
		ax[2].set_ylabel('Log Power', {'fontsize': ax_fs})

		# Adjust subplot spacing
		plt.subplots_adjust(hspace=0.4)

		print('Center-BW Corr: ' + str(corr_cen_bw) + ' with p of : ' + str(pcorr_cen_bw))
		print('Center-Power Corr: ' + str(corr_cen_pow) + ' with p of : ' + str(pcorr_cen_pow))
		print('BW-Power Corr: ' + str(corr_bw_pow) + ' with p of : ' + str(pcorr_bw_pow))


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

		# Set title for plots
		self.title = 'Group' 

	def add_subject(self, new_subj, add_all_oscs=True, add_vertex_bands=True):
		"""
		Input:
			new_subj	- MEG subject (instance of MegData)

		"""

		# Check if subject has data
		if(not new_subj.has_data):
			print("Empty meg data object. Cannot add data.")

		# Add All-Osc Data
		if(add_all_oscs):
			self.centers_all = np.append(self.centers_all, new_subj.centers_all)
			self.bws_all = np.append(self.bws_all, new_subj.bws_all)
			self.powers_all = np.append(self.powers_all, new_subj.powers_all)
			self.chis = np.append(self.chis, new_subj.chis)

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

		# Update subj count and subject number list
		self.nsubjs += 1
		self.subjs = np.append(self.subjs, new_subj.subnum)

		self.all_osc = add_all_oscs
		self.bands_vertex = add_vertex_bands


class MegProb():
	"""
	"""

	def __init__(self):
		pass


########################################################################
####################### OMEGAMAPPING - FUNCTIONS #######################
########################################################################


def clean_file_list(files_in, string):
    '''Takes a list of file names (strings), returns only those with 'string' in them.
    '''
    
    files_out = []

    for i in range(0, len(files_in)):
        if(string in files_in[i]):
            files_out.append(files_in[i])
            
    return files_out

def get_sub_nums(files_in):
	""" Takes a list of files. Returns a list of subject numbers. 
	"""

	sub_nums = []

	for f in files_in:
		str_split = f.split('_', 1)
		sub_nums.append(int(str_split[0]))

	return sub_nums

def get_cur_subj(cur_subj, files):
    '''Takes an int, and a list of file names (strings), returns the file name with the number in it.
    '''
    
    cur_subj_str = str(cur_subj)
    
    for i in range(0, len(files)):
        if(cur_subj_str in files[i]):
            return files[i]

def run_par_foof():
	""" NOT YET IMPLEMENTED.
	"""
	pass

def _get_osc(centers, powers, bws, osc_low, osc_high):
		""" Searches for an oscillations of specified frequency band. 
		Helper function for osc_per_vertex in MegData.
		"""

		osc_inds = (centers > osc_low) & (centers < osc_high)

		osc_cens = centers[osc_inds]
		osc_pows = powers[osc_inds]
		osc_bws = bws[osc_inds]

		cen, n = _get_single_osc(osc_cens)
		power, x = _get_single_osc(osc_pows)
		bw, x = _get_single_osc(osc_bws)

		return np.array([cen, power, bw, n])

def _get_single_osc(osc_in):
		'''Takes a vector, returns the first element, and the length.
		Helper function for osc_per_vertex in MegData.
		'''

		if(len(osc_in) == 0):
			return 0, 0
		elif(len(osc_in) == 1):
			return osc_in, 1
		else:
			return osc_in[0], len(osc_in)
