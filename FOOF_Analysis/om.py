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

		self.theta_low = 3
		self.theta_high = 8
		
		self.alpha_low = 8
		self.alpha_high = 13
		
		self.beta_low = 13
		self.beta_high = 30

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

		# Set boolean for whether contains data
		self.has_data = False
		self.anat_format = True

	def import_foof(self, subnum):
		""" Import FOOF results from pickle file. 
		"""

		# Set subject number for current data object
		self.subnum = subnum

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

	def osc_per_vertex(self, osc):
		"""
		"""

		## Initialize vectors to save results
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

		# Check if data in anatomy format
		if(not self.anat_format):
			print("Can't convert - subject not in anatomy format.")
			return
		
		# Set chis as empty so not accidentally used when data in this format
		self.chis = np.array([])

		# Flatten osc data into vectors. Uses C-style row-major order
		self.centers = self.centers.flatten('C')
		self.powers = self.powers.flatten('C')
		self.bws = self.bws.flatten('C')

		# Flattened vectors will have lots of zeros. Get only non-zero indices.
		non_zeros = np.nonzero(self.centers)
		self.centers = self.centers[non_zeros]
		self.powers = self.powers[non_zeros]
		self.bws = self.bws[non_zeros]

		# Get the number of oscillations
		self.nOscs = len(self.centers)

		# Update format
		self.anat_format = False

	def plot_all_oscs(self):
		""" Plots histogram distributions of oscillation centers, powers and bws. 
		"""

		# Check if in the right format
		if(self.anat_format):
			print("Can't plot - subject data in wrong format.")
			return

		# Plot Settings
		nBins = 100 	# Number of bins for histograms
		st_fs = 20 		# Super Title Font Size
		sp_fs = 18		# Subplit Title Font Size
		ax_fs = 16		# Axis Label Font Size

		# Set up subplots
		f, ax = plt.subplots(3, 1, figsize=(15, 15))
		# Set plot super-title
		plt.suptitle('Distributions of Oscillatory Parameters', fontsize=st_fs, fontweight='bold')

		# Subplot 1 - Center Frequency
		ax[0].hist(self.centers, nBins)
		ax[0].set_title('Center Frequency', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[0].set_xlabel('Frequency', {'fontsize': ax_fs})
		ax[0].set_ylabel('Count', {'fontsize': ax_fs})

		# Subplot 2 - Power
		ax[1].hist(np.log10(self.powers), nBins)
		ax[1].set_title('Oscillatory Power', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[1].set_xlabel('Log Power', {'fontsize': ax_fs})
		ax[1].set_ylabel('Count', {'fontsize': ax_fs})

		# Subplot 3 - Bandwidth
		ax[2].hist(self.bws, nBins)
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

		# Check if in the right format
		if(not self.anat_format):
			print("Can't plot - subject data in wrong format.")
			return

		# Plot Settings
		nBins = 100
		t_fs = 20
		ax_fs = 16

		plt.hist(self.chis, nBins)
		plt.title('Slopes', {'fontsize': t_fs, 'fontweight': 'bold'})
		plt.xlabel('Chi Parameter', {'fontsize': ax_fs, 'fontweight': 'bold'})
		plt.ylabel('Count', {'fontsize': ax_fs, 'fontweight': 'bold'})

	def plot_comparison(self):

		# Check if in the right format
		if(self.anat_format):
			print("Can't plot - subject data in wrong format.")
			return

		# Plot Settings
		st_fs = 20 		# Super Title Font Size
		sp_fs = 18		# Subplit Title Font Size
		ax_fs = 16		# Axis Label Font Size
		vis_alp = 0.1	# Opacity for plotting

		# Set up subplots
		f, ax = plt.subplots(3, 1, figsize=(15, 15))
		# Set plot super-title
		plt.suptitle('Oscillation Parameter Comparisons (Single Subject)', fontsize=st_fs, fontweight='bold')

		## Centers vs. Bandwidth
		# Check Correlation
		corr_cen_bw, pcorr_cen_bw = pearsonr(self.centers, np.log10(self.bws))

		# Plot
		ax[0].plot(self.centers, np.log10(self.bws), '.', alpha=vis_alp)
		ax[0].set_title('Center vs. Bandwidth', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[0].set_xlabel('Centers', {'fontsize': ax_fs})
		ax[0].set_ylabel('BW', {'fontsize': ax_fs})


		## Centers vs. Power
		# Check Correlation
		corr_cen_pow, pcorr_cen_pow = pearsonr(self.centers, np.log10(self.powers))

		# Plot
		ax[1].plot(self.centers, np.log10(self.powers), '.', alpha=vis_alp)
		ax[1].set_title('Center vs. Power', {'fontsize': sp_fs, 'fontweight': 'bold'})
		ax[1].set_xlabel('Centers', {'fontsize': ax_fs})
		ax[1].set_ylabel('Log Power', {'fontsize': ax_fs})


		## Bandwidth vs. Power
		# Check Correlation
		corr_bw_pow, pcorr_bw_pow = pearsonr(np.log10(self.bws), np.log10(self.powers))

		# Plot
		ax[2].plot(np.log10(self.bws), np.log10(self.powers), '.', alpha=vis_alp)		
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

		self.nsubjs = int()
		self.subjs = []

		self.chis = np.array([])
		self.centers = np.array([])
		self.powers = np.array([])
		self.bws = np.array([])

	def add_subject(self, new_subj):
		"""
		Input:
			new_subj	- MEG subject (instance of MegData)

		NOTE: THIS WON'T WORK RIGHT NOW, SHAPE OF RESULTS
		"""

		# Check if subject has data
		if(not new_subj.has_data):
			print("Empty meg data object. Cannot add.")

		# Add Data
		self.centers = np.append(self.centers, new_subj.centers)
		self.bws = np.append(self.bws, new_subj.bws)
		self.powers = np.append(self.powers, new_subj.powers)
		self.chis = np.append(self.chis, new_subj.chis)

		# Update subj count and subject number list
		self.nsubjs += 1
		self.subjs = np.append(self.subjs, new_subj.subnum)


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
