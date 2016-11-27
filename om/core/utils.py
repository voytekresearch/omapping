"""DOCSTRING"""

from __future__ import print_function, division

import numpy as np

from om.core.errors import InconsistentDataError

##
##
##

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
    files_out : list of str
        A list of the file names that contain the given string.
    """

    # Initialize variable of files to return
    files_out = []

    # Loop through given files, keeping those that contain string
    for i in range(len(files_in)):
        if string.lower() in files_in[i].lower():
            files_out.append(files_in[i])

    # Check if list is empty
    if not files_out:
        print('No files found!')

    return files_out


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
        Subject numbers for all files.
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

    Returns
    -------
    ind_st_x : int
        Starting index for the x axis data.
    ind_en_x : int
        Ending index for the x axis data.
    ind_st_y : int
        Starting index for the y axis data.
    ind_en_y : int
        Ending index for the y axis data.

    Raises
    ------
    InconsistentDataError
        xx

    Notes
    -----
    - For 'all', 'left', 'right', the 'a' and 'b' section indices are the same.
        The 'b' indices are different when getting the 'lr' and 'rl' comparisons.
    """

    # Set section indices
    if section is 'all':
        ind_st_x = ind_st_y = 0
        ind_en_x = ind_en_y = n_ROIs
    elif section is 'left':
        ind_st_x = ind_st_y = roi_lr.index('L')
        ind_en_x = ind_en_y = _find_last(roi_lr, 'L')
    elif section is 'right':
        ind_st_x = ind_st_y = roi_lr.index('R')
        ind_en_x = ind_en_y = _find_last(roi_lr, 'R')
    elif section is 'lr':
        ind_st_x = roi_lr.index('L')
        ind_en_x = _find_last(roi_lr, 'L')
        ind_st_y = roi_lr.index('R')
        ind_en_y = _find_last(roi_lr, 'R')
    elif section is 'rl':
        ind_st_x = roi_lr.index('R')
        ind_en_x = _find_last(roi_lr, 'R')
        ind_st_y = roi_lr.index('L')
        ind_en_y = _find_last(roi_lr, 'L')
    else:
        raise InconsistentDataError('Section range is unclear.')

    return ind_st_x, ind_en_x, ind_st_y, ind_en_y


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

def extract_foof_pickle(results):
    """Pull out the foof data from list.

    Parameters
    ----------
    results : list of tuple
        FOOF results - (slope (float), centers (1d array), amps (1d array), bws (1d array)).

    Returns
    -------
    centers : 2d array
        Matrix of all centers for all PSDs.
    powers : 2d array
        Matrix of all powers for all PSDs.
    bws : 2d array
        Matrix of all bws for all PSDs.
    slopes : 1d array
        Slope value for each PSD.
    n_psds : int
        The number of PSDs.
    """

    # Check how many psds there are
    n_psds = len(results)

    # Initialize numpy arrays to pull out different result parameters
    slopes = np.zeros([n_psds])
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

##
##
##

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
