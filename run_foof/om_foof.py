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

###############################
###### OM_FOOF Functions ######
###############################

def clean_file_list(files_in, string):
    """"   """
    
    files_out = []

    for i in range(0, len(files_in)):
        if(string in files_in[i]):
            files_out.append(files_in[i])
            
    return files_out


def load_meg_dat(meg_path, subj_num):
    """   """

    #
    mat_file = 'psd_source_median_' + str(subj_num)
    file_name = os.path.join((meg_path + 'Subject_' + str(subj_num)), mat_file)

    #
    data_mat = sio.loadmat(file_name, appendmat=True, struct_as_record=False, squeeze_me=True)

    # Pull out data from dictionary
    freqs = data_mat['Freqs']
    psd = data_mat['TF']

    # Label data is also available: Scout names if it's scout data, or Sensor/vertex numbers. 
    #labels = data_mat['RowNames']

    #
    return psd, freqs


def extract_psd(psd, freqs, f_low, f_high):
    """   """

    # Drop frequencies below f_low
    f_low_mask = freqs > f_low
    freqs_ext = freqs[f_low_mask]
    psd_ext = psd[:, f_low_mask]

    # Drop frequencies above f_high
    f_high_mask = freqs_ext < f_high
    freqs_ext = freqs_ext[f_high_mask]
    psd_ext = psd_ext[:, f_high_mask]

    return psd_ext, freqs_ext

def meg_foof(psd_ext, freqs_ext, method, min_p, freq_res):
    """   """

    # Check how many PSDs there are
    [nPSDs, nFreqs] = np.shape(psd_ext)

    # Initialize foof
    foof = FOOF(min_p=min_p, res=freq_res, fmin=freqs_ext.min(), fmax=freqs_ext.max())

    # Set up PSD as a list of 2-D np arrays
    psd_list = list(psd_ext)
    for i in range(0, nPSDs):
        psd_list[i] = np.reshape(psd_list[i], [len(freqs_ext), 1])

    #
    if method is 'linear':
        results = [_run_foof_l(foof, freqs_ext, psd) for psd in psd_list]

    #
    if method is 'parallel':
        pool = Pool(4)
        results = pool.map(_run_foof_p, psd_list)
        pool.close()
        pool.join()

    return results


def save_pickle(results, save_path, subj):
    """   """

    save_name = str(subj) + '_Foof_vertex.p'
    foof_save_path = os.path.join(save_path, save_name)
    pickle.dump(results, open(foof_save_path, 'wb'))

#######################################################
############## OM_FOOF Private Functions ##############
#######################################################

def _run_foof_l(foof, freqs_ext, psd_ext):
    """   """

    # Fit FOOF
    foof.model(freqs_ext, psd_ext)
    
    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)


def _run_foof_p(psd_ext):
    """   """

    # Fit FOOF
    foof.model(freqs_ext, psd_ext)
    
    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)




