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

#######################################################################
########################## OM_FOOF Functions ##########################
#######################################################################

def clean_file_list(files_in, string):
    """"Takes a list of files and returns only a specified set of files. 

    Parameters
    ----------
    files_in : list (str)
        A list of strings, each one being a file name.
    string : str
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
        xx
    freqs_ext : 1d array
        xx
    min_p : float
        xx
    freqs_res : float
        xx
    method : str
        Which method to use to run FOOF. 
        Options:
            'linear', 'parallel'

    Returns
    -------
    results : ?
        xx

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
        pool = Pool(4)
        results = pool.map(_run_foof_p, psd_list)
        pool.close()
        pool.join()

    return results


def save_pickle(results, save_path, sub_num):
    """Save out the FOOF results as a pickle file. 


    Parameters
    ----------
    results : ?
        xx
    save_path: str
        xx
    sub_num : int
        xx

    """

    # Set save name and path
    save_name = str(sub_num) + '_Foof_vertex.p'
    foof_save_path = os.path.join(save_path, save_name)
    
    # Save out data to pickle file
    pickle.dump(results, open(foof_save_path, 'wb'))


def save_csv(results, save_path, subj):
    """Save out the FOOF results as a csv file. 
    NOTE: Not yet implemented. 


    Parameters
    ----------

    Returns
    -------

    """

    pass

def load_pickle(save_path, sub_num):
    """   """
    
    subj_path = os.path.join()

    results = pickle.load(open(subj_path, 'rb'))
    return results

def load_csv():
    """   
    NOTE: Not yet implemented. 
    """
    pass

def conv_pickle_csv():
    """   
    NOTE: Not yet implemented. 
    """
    pass

    # Load pickle file

    # Convert to format for csv (??)

    # Save as csv file

def conv_csv_pickle():
    """
    NOTE: Not yet implemented. 
    """
    pass

    # Load csv file

    # Convert to format for pickle (??)

    # Save as pickle file


def get_sub_nums(files_in):
    """Takes a list of files. Returns a list of subject numbers.

    Parameters
    ----------
    files_in : list (str)
        List of filenames.
    
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
        subnums.append(int(str_split[0]))

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


################################################################
################## OM GEN - Private Functions ##################
################################################################

def _run_foof_l(foof, freqs_ext, psd_ext):
    """
    """

    # Fit FOOF
    foof.model(freqs_ext, psd_ext)
    
    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)


def _run_foof_p(psd_ext):
    """   
    """

    # Fit FOOF
    foof.model(freqs_ext, psd_ext)
    
    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)

