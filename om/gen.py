from __future__ import print_function, division
import os
import sys
import pickle
import numpy as np
import scipy.io as sio

# Import FOOF (use sys to add location to path, then import)
sys.path.append('/Users/thomasdonoghue/Documents/GitCode/')
from foof.fit import FOOF

###################################################################################
########################## OMEGAMAPPIN - GENERAL CLASSES ##########################
###################################################################################

class OMDB():
    """Class to hold database information for MEG project. """

    def __init__(self, dat_source='both'):

        # Check dat_source is acceptable
        pos_sources = ['both', 'OMEGA', 'HCP']
        if dat_source not in pos_sources:
            raise UnknownDataSourceError('Source to load data not understood.')

        # Save to object which data source is being used
        self.dat_source = dat_source

        # Set base path for OMEGA data
        self.project_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/'

        # Set paths for different data types
        self.maps_path = os.path.join(self.project_path, 'Maps')
        self.meg_path = os.path.join(self.project_path, 'MEG')
        self.corrs_path = os.path.join(self.project_path, 'Corrs')

        # Set paths for MEG data types
        self.psd_base_path = os.path.join(self.meg_path, '2-PSDs')
        self.foof_base_path = os.path.join(self.meg_path, '3-FOOF')
        self.viz_path = os.path.join(self.meg_path, '4-Viz')

        # Set path for database specific stuff
        if dat_source is 'both':
            self.psd_path = ''
            self.foof_path = ''
        else:
            self.psd_path = os.path.join(self.psd_base_path, dat_source)
            self.foof_path = os.path.join(self.foof_base_path, dat_source)

        # Set paths to save data out to
        self.processed_path = os.path.join(self.meg_path, '5-Processed')
        self.md_save_path = os.path.join(self.processed_path, 'md_pickle')
        self.mc_save_path = os.path.join(self.processed_path, 'mc_pickle')


    def check_dat_files(self, dat_type, save_type='pickle', verbose=True):
        """Checks what data files are available.

        Parameters
        ----------
        dat_type : str
            Which data type to check files for.
                Options: {'PSD', 'foof'}
        save_type : str, optional
            Which file type to check files for. Only used for foof files.
                Options: {'pickle', 'csv'}

        Returns
        -------
        sub_nums : list(int)
            A list of subject numbers of all the available files.
        """

        # Set up which files to look for
        if dat_type is 'PSD':
            dat_path = self.psd_base_path
            word = 'subject_'
            save_type = ''
            f_l = 'last'
        elif dat_type is 'foof':
            dat_path = self.foof_base_path
            word = 'foof'
            f_l = 'first'

        # ?
        if self.dat_source is not 'both':
            sub_nums = _check_files(os.path.join(dat_path, self.dat_source, save_type), word, f_l)
            source = [self.dat_source] * len(sub_nums)

        # ?
        else:
            sub_nums_omega = _check_files(os.path.join(dat_path, 'OMEGA', save_type), word, f_l)
            n_omega = len(sub_nums_omega)

            sub_nums_hcp = _check_files(os.path.join(dat_path, 'HCP', save_type), word, f_l)
            n_hcp = len(sub_nums_hcp)

            sub_nums = sub_nums_omega + sub_nums_hcp
            source = (['OMEGA'] * n_omega) + (['HCP'] * n_hcp)

        # If requested, print out the list of subject numbers
        if verbose:
            print('\nNumber of Subjects available: ' + str(len(sub_nums)) + '\n')
            print('Subject numbers with FOOF data available: \n' + str(sub_nums) + '\n')

        return sub_nums, source

    def check_res_files(self, res_type, verbose=True):
        """Checks what result files are available.

        Parameters
        ----------
        res_type : str
            Which data type to check files for.
                Options: {'md', 'mc'}

        Returns
        -------
        files : list(str)
            A list of all the available files.
        """

        # Settings
        word = 'Res'
        
        # Set up which files to look for
        if res_type is 'md':
            dat_path = self.md_save_path
        elif res_type is 'mc':
            dat_path = self.mc_save_path

        # Get files
        files = os.listdir(dat_path)
        files = clean_file_list(files, word)

        # If requested, print out the list of subject numbers
        if verbose:
            print('\nNumber of files available: ' + str(len(files)) + '\n')
            print('Files available: \n' + ('\n'.join(files)) + '\n')

        # Return file list
        return files


class Osc():
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

"""
TODO: Update Osc and everything that calls it to use a more flexible dictionary mapping. 
This is so it's much easier to change / add / remove oscillation bands. 

Something like: 

class Osc():

    def __init__(self, input_bands=None):

        # If supplied, use the given dictionary of oscillation bands. 
        if input_bands:
            self.bands = input_bands

        # Otherwise, use the default oscillation bands.
        else:
            self.bands = dict({'Theta': (3, 8), 
                               'Alpha': (8, 13), 
                               'Beta': (13, 30), 
                               'LowGamma': (30, 40)})
"""

class FigInfo():
    """Object to hold settings to save figures. """

    def __init__(self):

        # Default Settings - font sizes
        self.t_fs = 22           # Title font size
        self.sp_fs = 20          # Subplot title font size
        self.ax_fs = 20          # Axis font size
        self.ti_fs = 14          # Ticks font size

        # Default Settings - other settings
        self.ax_lw = 2.5

        # Default Settings - what to add to plot
        self.add_title = False

        # Plot Information
        self.title = 'Group'
        self.vis_opac = 0.005

        # Save Information
        self.save_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/4-Figures/MegAnalysis/'
        self.format = 'svg'
        self.bbox = 'tight'
        self.dpi = 150


###################################################################################
################################ OM_FOOF Functions ################################
###################################################################################

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
        if string.lower() in files_in[i].lower():
            files_out.append(files_in[i])

    # Check if list is empty
    if not files_out:
        print('No files found!')

    return files_out


def load_meg_psds(dat_source, meg_path, subj_num):
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
    if dat_source is 'OMEGA':
        mat_file = 'psd_source_median_' + str(subj_num)
    elif dat_source is 'HCP':
        mat_file = 'PSD_Source_Median_' + str(subj_num)
    file_name = os.path.join(meg_path, ('Subject_' + str(subj_num)), mat_file)

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
    NOTE: OLD - NOT CURRENTLY USED.
        REPLACED BY STANDALONE SCRIPT.

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
            osc_csv.write((", ".join(repr(el) for el in cur_osc_dat)) + '\n')


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
    files_in : list[str]
        List of filenames.
    f_l : str
        Whether subject numbers are first or last in file name.
        Options: 'first', 'last'

    Returns
    -------
    subnums : list[int]
        List of subject numbers.
    """

    # Check and remove file extensions
    files_in = rm_files_ext(files_in)

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
    files : list[str]
        List of files to search through.

    Returns
    -------
    subj_file : str
        File name of specific subject's file.
    """

    # Make sure given number is a string
    cur_subj_str = str(cur_subj)

    # Loop through files, looking for one which contains #
    for i in range(0, len(files)):
        if cur_subj_str in files[i]:
            return files[i]


def rm_files_ext(files_in):
    """Removes file extensions for list of files given.

    Parameters
    ----------
    files_in : list[str]
        A list of file and/or directory names.

    Returns
    -------
    files_out : list[str]
        A list of file and/or directory names with file extensions removed.
    """

    # Initialize list to store output
    files_out = []

    # Loop through all files
    for f_name in files_in:

        # If there is an extension, remove it
        if '.' in f_name:
            cur_file = f_name.split('.')[0]
        else:
            cur_file = f_name

        # Gather all file names to return
        files_out.append(cur_file)

    return files_out


def get_section(section, nROIs, roi_lr):
    """Get indices for desired section of ...

    Parameters
    ----------
    section : str
        Which section to get indices for.
    nROIs : int
        The number of ROIs.
    roi_lr : list(str)
        List of L/R for each ROI.
    """

    # Set section indices
    if section is 'all':
        ind_st_a = ind_st_b = 0
        ind_en_a = ind_en_b = nROIs
    elif section is 'left':
        ind_st_a = ind_st_b = roi_lr.index('L')
        ind_en_a = ind_en_b = _find_last(roi_lr, 'L')
    elif section is 'right':
        ind_st_a = ind_st_b = roi_lr.index('R')
        ind_en_a = ind_en_b =  _find_last(roi_lr, 'R')
    elif section is 'lr':
        ind_st_a = roi_lr.index('L')
        ind_en_a = _find_last(roi_lr, 'L')
        ind_st_b = roi_lr.index('R')
        ind_en_b = _find_last(roi_lr, 'R')
    elif section is 'rl':
        ind_st_a = roi_lr.index('R')
        ind_en_a = _find_last(roi_lr, 'R')
        ind_st_b = roi_lr.index('L')
        ind_en_b = _find_last(roi_lr, 'L')
    else:
        print('Section range unclear!')
        ind_st_a = ind_en_a = ind_st_b = ind_en_b = 0

    # Return indices
    return ind_st_a, ind_en_a, ind_st_b, ind_en_b


#########################################################################################
###################################### OM GEN - ERRORS ##################################
#########################################################################################

class UnknownDataSourceError(Exception):
    """An Error indicating data source specification is not understood."""
    pass

class UnknownDataTypeError(Exception):
    """An Error indicating data type specification is not understood."""
    pass

########################################################################################
############################## OM GEN - Private Functions ##############################
########################################################################################

def _run_foof_l(foof, freqs_ext, psd_ext):
    """Local helper function to run FOOF linearly.

    Used by meg_foof().
    NOTE: CURRENTLY UNUSED.

    Parameters
    ----------
    foof : FOOF() object
        FOOF object to model 1/f & oscillations.
    freqs_ext : 1d array
        Vector of frequency values for the psd.
    psd_ext : 1d array
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
    NOTE: CURRENTLY UNUSED.
    """

    # Fit FOOF
    foof.model(freqs_ext, psd_ext)

    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)


def _check_files(path, word, f_l):
    """???

    Parameters
    ----------
    path : ?
        xx
    word : ?
        xx
    f_l : ?
        xx

    Returns
    -------
    sub_nums : ?
        xx
    """

    files = os.listdir(path)
    files = clean_file_list(files, word)

    sub_nums = get_sub_nums(files, f_l)

    return sub_nums


def _find_last(input_list, wanted):
    """Find the index of the last instance of a wanted element.

    Parameters
    ----------
    input_list : list
        A list to search through.
    wanted : str
        The element for which the last index is wanted.

    From here: http://stackoverflow.com/questions/6890170/how-to-find-the-last-occurrence-of-an-item-in-a-python-list
    """

    # Run through the list, backwards
    for ind, el in enumerate(reversed(input_list)):

        # If element is the wanted one, return index
        if el == wanted:
            return len(input_list) - 1 - ind
