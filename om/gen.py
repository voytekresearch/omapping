from __future__ import print_function, division
import os
import sys
import pickle
import numpy as np
import scipy.io as sio

# Import FOOF (use sys to add location to path, then import)
sys.path.append('/Users/thomasdonoghue/Documents/GitCode/')
from foof.fit import FOOF

#####################################################################################
########################## OMEGAMAPPIN - GENERAL - CLASSES ##########################
#####################################################################################

class OMDB(object):
    """Class to hold database information for MEG project. """

    def __init__(self, dat_source='both'):
        """

        Parameters
        ----------
        dat_source : {'both', 'OMEGA', 'HCP'}
            xx
        """

        # Check dat_source is acceptable
        pos_sources = ['both', 'OMEGA', 'HCP']
        if dat_source not in pos_sources:
            raise UnknownDataSourceError('Source to load data not understood.')

        # Save to object which data source is being used
        self.dat_source = dat_source

        # Set base path for OMEGA data
        self.project_path = ("/Users/thomasdonoghue/Documents/Research/"
                             "1-Projects/OMEGA/2-Data/")

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
        dat_type : {'PSD', 'foof'}
            Which data type to check files for.
        save_type : {'pickle', 'csv'}, optional (default = 'pickle')
            Which file type to check files for. Only used for foof files.
        verbose : boolean, optional
            Whether to print out information during run.

        Returns
        -------
        sub_nums : list of int
            A list of subject numbers of all the available files.
        source : list of str
            A list with the source database for each available file.
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

        # If looking for a particular database, find file, get subject numbers and source
        if self.dat_source is not 'both':
            sub_nums = _check_files(os.path.join(dat_path, self.dat_source, save_type), word, f_l)
            source = [self.dat_source] * len(sub_nums)

        # If looking across both databases, get info from each database and then combine
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
        res_type : {'md', 'mc'}
            Which data type to check files for.
        verbose : boolean, optional (default = True)
            Whether to print out information during run.

        Returns
        -------
        files : list of str
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


class Osc(object):
    """Class to hold definition of oscillation bands."""

    def __init__(self, default=False, input_bands=None):
        """Initialize the Osc() object.

        Parameters
        ----------
        default : boolean, optional (default = False)
            Whether to use the default oscillation bands. Default it False.
        input_bands : dict, optional (default = None)
            A dictionary of oscillation bands to use.

        Notes:
        - If supplied, an input_bands will over-ride the default bands option,
            even if it is set as True.
        """

        # Initialize bands as a dictionary
        self.bands = dict()

        # If requested use the default oscillation bands
        if default:
            self.bands = dict({'Theta': (3, 8),
                               'Alpha': (8, 13),
                               'Beta': (13, 30),
                               'LowGamma': (30, 40)})

        # If supplied, use the given dictionary of oscillation bands
        if input_bands:
            self.bands = input_bands


    def add_band(self, band_name, band_lims):
        """Add a new oscillation band definition.

        Parameters
        ----------
        band_name : str
            The name of the new oscillation band.
        band_lims : tuple(float, float)
            The lower and upper frequency limit of the band.
        """

        # Safety check that limits are in correct order
        if not band_lims[0] < band_lims[1]:
            raise InconsistentDataError('Band limits are incorrect.')

        # Add the given band to oscillation definition
        self.bands[band_name] = band_lims


    def rm_band(self, old_band):
        """Remove a previously defined oscillation band.

        Parameters
        ----------
        old_band : str
            Band name to remove from oscillation band definitions.
        """

        # Remove requested band from oscillation definition
        self.bands.pop(old_band)


class FigInfo(object):
    """Object to hold settings to save figures. """

    def __init__(self):

        # Default Settings - font sizes
        self.t_fs = 22         # Title font size
        self.sp_fs = 20        # Subplot title font size
        self.ax_fs = 20        # Axis font size
        self.ti_fs = 14        # Ticks font size

        # Default Settings - other settings
        self.ax_lw = 2.5

        # Default Settings - what to add to plot
        self.add_title = False

        # Plot Information
        self.title = 'Group'
        self.vis_opac = 0.005

        # Save Information
        self.save_path = ("/Users/thomasdonoghue/Documents/Research/"
                          "1-Projects/OMEGA/4-Figures/MegAnalysis/")
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
    files_in : list of str
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
    dat_source : {'OMEGA', 'HCP'}
        xx
    meg_path : str
        Path where data to load is located.
    subj_num : int
        Subject identifier number.

    Returns
    -------
    psd : 2d array
        Matrix of PSD results for each vertex [n_verts, n_freqs].
    freqs : 1d array
        Vector of the frequencies estimated in the PSD.
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

    return psd, freqs


def extract_psd(psd, freqs, f_low, f_high):
    """Extract frequency range of interest from PSD data.

    Parameters
    ----------
    psd : 2d array
        Matrix of PSD results for each vertex [n_verts, n_freqs].
    freqs : 1d array
        Vector of the frequencies estimated in the PSD.
    f_low : float
        Lower bound of frequencies to extract.
    f_high : float
        Upper bound of frequencies to extract.

    Returns
    -------
    psd_ext : 2d array
        Matrix of extracted PSD results for each vertex [n_verts, n_freqs].
    freqs_ext : 1d array
        Vector of extracted frequencies estimated in the PSD.
    """

    # Drop frequencies below f_low
    f_low_mask = freqs >= f_low
    freqs_ext = freqs[f_low_mask]
    psd_ext = psd[:, f_low_mask]

    # Drop frequencies above f_high
    f_high_mask = freqs_ext <= f_high
    freqs_ext = freqs_ext[f_high_mask]
    psd_ext = psd_ext[:, f_high_mask]

    return psd_ext, freqs_ext


def save_foof_pickle(results, save_path, sub_num):
    """Save out the FOOF results as a pickle file.

    Parameters
    ----------
    results : list of tuple
        List of tuples of FOOF results. Length of list is number of vertices.
            Tuple is: (slope (float), centers (1d array), amps (1d array), bws (1d array)).
    save_path: str
        Filepath of where to save out the file.
    sub_num : int
        Subject identifier number.
    """

    # Set save name and path
    save_name = str(sub_num) + '_Foof_Vertex.p'
    foof_save_path = os.path.join(save_path, 'pickle', save_name)

    # Save out data to pickle file
    pickle.dump(results, open(foof_save_path, 'wb'))


def save_foof_csv(results, save_path, sub_num):
    """Save out the FOOF results as a csv file.

    Parameters
    ----------
    results : list of tuple
        List of tuples of FOOF results. Length of list is number of vertices.
            Tuple is: (slope (float), centers (1d array), amps (1d array), bws (1d array)).
    save_path : str
        Filepath of where to save out the file.
    sub_num : int
        Subject identifier number.
    """

    # Set index values for oscillation parameters
    i_cen = 1
    i_amp = 2
    i_bw = 3

    # Check how many vertices there are
    n_verts = len(results)

    # Initialize file names for oscillation and slope csv files
    csv_sl_fname = save_path + '/csv/' + str(sub_num) + '_Slopes.csv'
    csv_osc_fname = save_path + '/csv/' + str(sub_num) + '_Oscs.csv'

    # Open files to write to
    sl_csv = open(csv_sl_fname, 'w')
    osc_csv = open(csv_osc_fname, 'w')

    # Loop through each vertex
    for vert in range(n_verts):

        # Save out slope value to csv file
        sl_csv.write(str(results[vert][0]) + '\n')

        # Check how oscillations at current vertex
        n_oscs = len(results[vert][1])

        # Loop through each oscillation, writing each one to file
        for osc in range(n_oscs):

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
    results : list of tuple
        List of tuples of FOOF results. Length of list is number of vertices.
            Tuple is: (slope (float), centers (1d array), amps (1d array), bws (1d array)).
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


def extract_foof_pickle(results):
    """Pull out the foof data from pickled list.

    Parameters
    ----------
    results : ?
        xx

    Returns
    -------
    centers : 2d array
        xx
    powers : 2d array
        xx
    bws : 2d array
        xx
    slopes : 1d array
        xx
    n_psds : int
        xx
    """

    # Check how many psds there are
    n_psds = len(results)

    # Initialize numpy arrays to pull out different result parameters
    slopes = np.zeros([n_psds, 1])
    centers = np.zeros([n_psds, 8])
    powers = np.zeros([n_psds, 8])
    bws = np.zeros([n_psds, 8])

    # Pull out the data from each vertex
    for i in range(n_psds):
        slopes[i] = results[i][0]
        centers[i, 0:len(results[i][1])] = results[i][1]
        powers[i, 0:len(results[i][2])] = results[i][2]
        bws[i, 0:len(results[i][3])] = results[i][3]

    return centers, powers, bws, slopes, n_psds


def load_foof_csv():
    """
    NOTE: Not yet implemented.
    """

    pass


def get_sub_nums(files_in, f_l):
    """Takes a list of files. Returns a list of subject numbers.

    Parameters
    ----------
    files_in : list of str
        List of filenames.
    f_l : {'first', 'last'}
        Whether subject numbers are first or last in file name.

    Returns
    -------
    subnums : list of int
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
    files : list of str
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
    files_in : list of str
        A list of file and/or directory names.

    Returns
    -------
    files_out : list of str
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


def get_section(section, n_ROIs, roi_lr):
    """Get indices for desired section of connectivity matrix.

    Parameters
    ----------
    section : {'all', 'left', 'right', 'lr', 'rl'}
        Which section to get indices for.
    n_ROIs : int
        The number of ROIs.
    roi_lr : list of str
        List of L/R for each ROI.
    """

    # Set section indices
    if section is 'all':
        ind_st_a = ind_st_b = 0
        ind_en_a = ind_en_b = n_ROIs
    elif section is 'left':
        ind_st_a = ind_st_b = roi_lr.index('L')
        ind_en_a = ind_en_b = _find_last(roi_lr, 'L')
    elif section is 'right':
        ind_st_a = ind_st_b = roi_lr.index('R')
        ind_en_a = ind_en_b = _find_last(roi_lr, 'R')
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
        raise InconsistentDataError('Section range is unclear.')

    # Return indices
    return ind_st_a, ind_en_a, ind_st_b, ind_en_b


def meg_foof(psd_ext, freqs_ext, min_p, freq_res):
    """Run FOOF on MEG-PSD data.

    NOTE: This is an old version, can run foof linearly.
        For parallel, use standalone script.

    Parameters
    ----------
    psd_ext : 2d array
        Matrix of PSDs in the form of [n_verts, n_freqs].
    freqs_ext : 1d array
        Vector of the frequency values for each power value in psd_ext.
    min_p : float
        Minimum probability for splitting peaks. Parameter for FOOF.
    freqs_res : float
        Frequency resolution.

    Returns
    -------
    results : list of tuple
        List of tuples of FOOF results. Length of list is number of vertices.
            Tuple is: (slope (float), centers (1d array), amps (1d array), bws (1d array)).
    """

    # Check how many PSDs there are
    [n_PSDs, n_freqs] = np.shape(psd_ext)

    # Initialize foof
    foof = FOOF(min_p=min_p, res=freq_res, fmin=freqs_ext.min(), fmax=freqs_ext.max())

    # Set up PSD as a list of 2-D np arrays
    psd_list = list(psd_ext)
    for i in range(0, n_PSDs):
        psd_list[i] = np.reshape(psd_list[i], [len(freqs_ext), 1])

    # Run FOOF linearly
    results = [_run_foof_l(foof, freqs_ext, psd) for psd in psd_list]

    return results


#########################################################################################
###################################### OM GEN - ERRORS ##################################
#########################################################################################

class UnknownDataSourceError(Exception):
    """An Error indicating data source specification is not understood."""
    pass

class UnknownDataTypeError(Exception):
    """An Error indicating data type specification is not understood."""
    pass

class InconsistentDataError(Exception):
    """An Error indicating there is a fatal inconsistency in data."""
    pass

class DataNotComputedError(Exception):
    """An Error indicating some required data has not been computed."""
    pass

########################################################################################
############################## OM GEN - Private Functions ##############################
########################################################################################

def _run_foof_l(foof, freqs_ext, psd_ext):
    """Local helper function to run FOOF linearly.

    Used by meg_foof().

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
            Tuple is (slope (float), centers (1d array), amps (1d array), bws (1d array)).
    """

    # Fit FOOF
    foof.model(freqs_ext, psd_ext)

    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)


def _check_files(path, word, f_l):
    """Checks a directory, getting desired files and returning subject numbers.

    Parameters
    ----------
    path : str
        Path to directory to examine.
    word : str
        Word to search for in file names to keep.
    f_l : {'first', 'last'}
        Whether subject number is at start or end of file name.

    Returns
    -------
    sub_nums : list of int
        A list of subject numbers of all the available files.
    """

    # Get list of files in desired directory
    files = os.listdir(path)
    files = clean_file_list(files, word)

    # Get the list of subject numbers from directory
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
    """

    # Run through the list, backwards
    for ind, elem in enumerate(reversed(input_list)):

        # If element is the wanted one, return index
        if elem == wanted:
            return len(input_list) - 1 - ind
