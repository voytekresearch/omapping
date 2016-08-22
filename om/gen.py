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

#
#
#

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


################################################################
################## OM GEN - Private Functions ##################
################################################################

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


##########################################################################################
###############################  OM CL - PRIVATE FUNCTIONS  ##############################
##########################################################################################


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

    #
    osc_inds = (centers > osc_low) & (centers < osc_high)
    osc_cens = centers[osc_inds]

    return osc_cens


def _get_demo_csv(subnum):
    """Get demographic information from csv file for specified subject.

    Parameters
    ----------
    subnum : int
        Subject number to get demographic info for. 

    Returns
    -------
    sex : str
        Sex ['M'/'F'] of specified subject. 
    age : int
        Age (in whole years) of specified subject. 
    """

    # Set up paths for demographic info csv file
    csv_data_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/'
    csv_file_name = '00-Collin_Subjects.csv'
    csv_file = os.path.join(csv_data_path, csv_file_name)

    # Open csv file, loop through looking for right row, grab age & sex information
    with open(csv_file, 'rb') as f_name:
        reader = csv.reader(f_name, delimiter=',')
        for row in reader:
            if row[1] == str(subnum):
                sex = row[4]
                age = int(row[7])
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
            Options: 'mean', 'median'

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


def _get_map_names(names_file, path):
    """

    Parameters
    ----------
    names_files : ??
        XX
    path : ?
        XX

    Returns
    -------
    names : ?
        XX
    """

    # Get path to csv file
    csv_path = os.path.join(path, names_file)

    # Open csv file
    with open(csv_path, 'rb') as f_name:

        reader = csv.reader(f_name, delimiter=',')

        # Get list of names from first row in csv
        names = list(reader)[0]

    return names


def _init_meg_map_dict(length=0):
    """   """
    if length == 0:
        meg_map = dict([('Theta',     np.array([])),
                        ('Alpha',     np.array([])),
                        ('Beta',      np.array([])),
                        ('LowGamma',  np.array([])),
                        ('Slopes',    np.array([]))
                        ])

    else:
        meg_map = dict([('Theta',     np.zeros(length)),
                        ('Alpha',     np.zeros(length)),
                        ('Beta',      np.zeros(length)),
                        ('LowGamma',  np.zeros(length)),
                        ('Slopes',    np.zeros(length))
                        ])

    return meg_map


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