from __future__ import print_function, division
import os
import pickle
import numpy as np
import scipy.io as sio

# Import FOOF (use sys to add location to path, then import)
import sys
sys.path.append('/Users/thomasdonoghue/Documents/GitCode/')
from foof import syn
from foof.fit import FOOF

# Import Parallelization Packages
from multiprocessing import Pool, freeze_support

###########################################################################
###################### OMEGAMAPPIN - GENERAL CLASSES ######################
###########################################################################

class Osc:
    """Class to store oscillations parameters."""
    
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


#######################################################################
########################## OM_FOOF Functions ##########################
#######################################################################

def clean_file_list(files_in, string):
    """"Takes a list of files and returns only a specified set of files. 

    Parameters
    ----------
    files_in : list (str)
        A list of strings, each one being a file name.
    string : str OR list
        A string to look for in file list, to keep those who have it.

    Returns
    -------
    files_out : list (str)
        A list of the file names that contain the given string.
    """
    
    # Initialize variable of files to return
    files_out = []

    # Loop through given files, keeping those that contain string
    for i in range(0, len(files_in)):
        if(string in files_in[i]):
            files_out.append(files_in[i])

    # Check if list is empty
    if not files_out:
        print('No files found!')
            
    return files_out


def load_meg_psds(meg_path, subj_num):
    """Loads a requested subject's PSD-MEG data.

    Parameters
    ----------
    meg_path : str
        Path where data to load is located. 
    subj_num : int
        Subject identifier number. 

    Returns
    -------
    psd : 2d array
        xx
    freqs : 1d array
        xx
    """

    # Set file name and get full path
    mat_file = 'psd_source_median_' + str(subj_num)
    file_name = os.path.join((meg_path + 'Subject_' + str(subj_num)), mat_file)

    # Load MEG PSD data from matfile
    data_mat = sio.loadmat(file_name, appendmat=True, struct_as_record=False, squeeze_me=True)

    # Pull out data from dictionary
    freqs = data_mat['Freqs']
    psd = data_mat['TF']

    # Label data is also available: Scout names if it's scout data, or Sensor/vertex numbers. 
    #labels = data_mat['RowNames']

    return psd, freqs


def extract_psd(psd, freqs, f_low, f_high):
    """Extract frequency range of interest from PSD data. 

    Parameters
    ----------
    psd : 2d array
        xx
    freqs : 1d array
        xx
    f_low : float
        Lower bound of frequencies to extract. 
    f_high : float
        Upper bound of frequencies to extract. 

    Returns
    -------
    psd_ext : 2d array
        xx
    freqs_ext : 1d array
        xx
    """

    # Drop frequencies below f_low
    f_low_mask = freqs > f_low
    freqs_ext = freqs[f_low_mask]
    psd_ext = psd[:, f_low_mask]

    # Drop frequencies above f_high
    f_high_mask = freqs_ext < f_high
    freqs_ext = freqs_ext[f_high_mask]
    psd_ext = psd_ext[:, f_high_mask]

    return psd_ext, freqs_ext

def meg_foof(psd_ext, freqs_ext, min_p, freq_res, method):
    """Run FOOF on MEG-PSD data. 

    Parameters
    ----------
    psd_ext : 2d array
        Matrix of PSDs in the form of [nVerts, nFreqs].
    freqs_ext : 1d array
        Vector of the frequency values for each power value in psd_ext. 
    min_p : float
        Minimum probability for splitting peaks. Parameter for FOOF. 
    freqs_res : float
        Frequency resolution. 
    method : str
        Which method to use to run FOOF. 
        Options:
            'linear', 'parallel'

    Returns
    -------
    results : list
        List of tuples of FOOF results. Length of list is number of vertices. 
            Each tuple is (slope value, centers (list), amps (list), bws (list)).
    """

    # Check how many PSDs there are
    [nPSDs, nFreqs] = np.shape(psd_ext)

    # Initialize foof
    foof = FOOF(min_p=min_p, res=freq_res, fmin=freqs_ext.min(), fmax=freqs_ext.max())

    # Set up PSD as a list of 2-D np arrays
    psd_list = list(psd_ext)
    for i in range(0, nPSDs):
        psd_list[i] = np.reshape(psd_list[i], [len(freqs_ext), 1])

    # Run FOOF linearly
    if method is 'linear':
        results = [_run_foof_l(foof, freqs_ext, psd) for psd in psd_list]

    # Run FOOF in parallel
    if method is 'parallel':

        # Set up Pool and run
        pool = Pool(4)
        results = pool.map(_run_foof_p, psd_list)
        #results = pool.map(_run_foof_l, [foof, psd_list, freqs_ext])
        pool.close()
        pool.join()

    return results


def save_foof_pickle(results, save_path, sub_num):
    """Save out the FOOF results as a pickle file. 

    Parameters
    ----------
    results : list
        xx
    save_path: str
        Filepath of where to save out the file. 
    sub_num : int
        Subject identifier number. 
    """

    # Set save name and path
    save_name = str(sub_num) + '_Foof_Vertex.p'
    foof_save_path = os.path.join(save_path, save_name)
    
    # Save out data to pickle file
    pickle.dump(results, open(foof_save_path, 'wb'))


def save_foof_csv(results, save_path, sub_num):
    """Save out the FOOF results as a csv file. 

    Parameters
    ----------
    results : list
        xx
    save_path : str
        xx
    sub_num : int
        Subject identifier number. 
    """

    #
    i_cen = 1
    i_amp = 2
    i_bw = 3

    #
    nVert = len(results)

    #
    csv_sl_fname = save_path + '/csv/' + str(sub_num) + '_Slopes.csv'
    csv_osc_fname = save_path + '/csv/' + str(sub_num) + '_Oscs.csv'

    #
    sl_csv = open(csv_sl_fname, 'w')
    osc_csv = open(csv_osc_fname, 'w')

    #
    for vert in range(0, nVert):

        #
        sl_csv.write(str(results[vert][0]) + '\n')

        #
        nOscs = len(results[vert][1])

        # 
        for osc in range(0, nOscs):

            cur_osc_dat = list([vert + 1, results[vert][i_cen][osc], 
                results[vert][i_amp][osc], results[vert][i_bw][osc]])
            osc_csv.write( (", ".join( repr(el) for el in cur_osc_dat )) + '\n')


def load_foof_pickle(dat_path, sub_num):
    """Load FOOF data from pickle file. 

    Parameters
    ----------
    dat_path : str
        File name for where data is stored to load from. 
    sub_num : int
        Subject identifier number. 

    Returns
    -------
    results : ?
        xx
    """
    
    # Get list of available files to load
    files = os.listdir(os.path.join(dat_path, 'pickle'))
    files = clean_file_list(files, 'Foof_Vertex')

    # Get specific subject file
    cur_subj_file = get_cur_subj(sub_num, files)
    subj_path = os.path.join(dat_path, 'pickle', cur_subj_file)

    # Load file
    results = pickle.load(open(subj_path, 'rb'))
    
    return results


def load_foof_csv():
    """   
    NOTE: Not yet implemented. 
    """

    pass


def get_sub_nums(files_in, f_l):
    """Takes a list of files. Returns a list of subject numbers.

    Parameters
    ----------
    files_in : list (str)
        List of filenames.
    f_l : str
        Whether subject numbers are first or last in file name. 
        Options: 'first', 'last'
    
    Returns
    -------
    subnums : list (int)
        List of subject numbers. 
    """

    # Intialize variable to store subject numbers
    subnums = []

    # Loop through files, extracting subject numbers
    for f_name in files_in:
        str_split = f_name.split('_', 1)
        if f_l is 'first':
            subnums.append(int(str_split[0]))
        elif f_l is 'last':
            subnums.append(int(str_split[1]))

    return subnums


def get_cur_subj(cur_subj, files):
    """Returns the file name with the given number in it.

    Parameters
    ----------
    cur_subj : int
        Subject number to search for in given file list.
    files : list (str)
        List of files to search through. 

    Returns
    -------
    subj_file
        File name of specific subject's file. 
    """

    # Make sure given number is a string
    cur_subj_str = str(cur_subj)

    # Loop through files, looking for one which contains #
    for i in range(0, len(files)):
        if cur_subj_str in files[i]:
            return files[i]


########################################################################
###################### OM GEN - Private Functions ######################
########################################################################

def _run_foof_l(foof, freqs_ext, psd_ext):
    """Local helper function to run FOOF linearly. 

    Used by meg_foof(). 

    Parameters
    ----------
    foof : FOOF() object
        FOOF object to model 1/f & oscillations. 
    freqs_ext : 1d vector
        Vector of frequency values for the psd. 
    psd_ext : 1d vector
        Vector of power values for the psd. 

    Returns
    -------
    out : tuple
        Tuple of FOOF results. 
            Tuple is (slope value, centers (list), amps (list), bws (list)).
    """

    # Fit FOOF
    foof.model(freqs_ext, psd_ext)
    
    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)


def _run_foof_p(psd_ext):
    """ ???
    """

    # Fit FOOF
    foof.model(freqs_ext, psd_ext)
    
    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)
