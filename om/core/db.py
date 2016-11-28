"""DOCSTRING"""
from __future__ import print_function, division

import os

from om.core.utils import clean_file_list, get_sub_nums

#####################################################################################
############################## OMEGAMAPPIN - CORE - DB ##############################
#####################################################################################

class OMDB(object):
    """Class to hold database information for MEG project.

    Attributes
    ----------
    internal_path : str
        xx
    external_path : str
        xx
    maps_path : str
        Path to Maps data.
    corrs_path : str
        Path to Corrs data.
    md_save_path : str
        Path to save md data.
    mc_save_path : str
        Path to save mc data.
    meg_path : str
        Path to MEG data.
    psd_path : str
        Path to PSD data.
    foof_path : str
        Path to FOOF data.
    viz_path : str
        Path to vizualization data.
    """

    def __init__(self, auto_gen=True):
        """   """

        # Initialize base paths
        self.internal_path = ("/Users/thomasdonoghue/Documents/Research/"
                              "1-Projects/OMEGA/2-Data/OMData/")
        self.external_path = ("/Users/thomasdonoghue/Documents/Research/"
                              "1-Projects/OMEGA/2-Data/ExternalData/")

        # Initialize all internal paths
        self.maps_path = str()
        self.corrs_path = str()
        self.md_save_path = str()
        self.mc_save_path = str()

        # Initialize all external paths
        self.meg_path = str()
        self.psd_path = str()
        self.foof_path = str()
        self.viz_path = str()

        # Generate project paths
        if auto_gen:
            self.gen_paths()


    def gen_paths(self):
        """Generate all the full paths for the OM project."""

        # Set up internal data paths for Maps & Corrs
        self.maps_path = os.path.join(self.internal_path, 'Maps')
        self.corrs_path = os.path.join(self.internal_path, 'Corrs')

        # Set up internal paths to save data out to
        processed_path = os.path.join(self.internal_path, 'Processed')
        self.md_save_path = os.path.join(processed_path, 'meg')
        self.mc_save_path = os.path.join(processed_path, 'maps')

        # Set up external data paths
        self.meg_path = os.path.join(self.external_path, 'MEG')
        self.psd_path = os.path.join(self.meg_path, 'PSDs')
        self.foof_path = os.path.join(self.meg_path, 'FOOF')
        self.viz_path = os.path.join(self.meg_path, 'Viz')


    def check_dat_files(self, dat_type, dat_source='both', save_type='pickle', verbose=True):
        """Checks what data files are available.

        Parameters
        ----------
        dat_type : {'PSD', 'foof'}
            Which data type to check files for.
        dat_source : {'OMEGA', 'HCP', 'both'}
            Which database to check files for.
        save_type : {'pickle', 'csv'}, optional (default = 'pickle')
            Which file type to check files for. Only used for foof files.
        verbose : boolean, optional (default = True)
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
            dat_path = self.psd_path
            word = 'subject_'
            save_type = ''
            f_l = 'last'
        elif dat_type is 'foof':
            dat_path = self.foof_path
            word = 'foof'
            f_l = 'first'

        # If looking for a particular database, find file, get subject numbers and source
        if dat_source is not 'both':
            sub_nums = _check_files(os.path.join(dat_path, dat_source, save_type), word, f_l)
            source = [dat_source] * len(sub_nums)

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

        return files


    def check_map_files(self, verbose=True, return_files=False):
        """Gets the list of files in the map directories. Can return and/or print.

        Parameters
        ----------
        verbose : boolean, optional (default = True)
            Whether or not to print out available file names.
        return_files : boolean, optional (default = False)
            Whether or not to return lists of filenames.

        Returns
        -------
        osc_files : list of str
            All available oscillation files.
        slope_files : list of str
            All available slope files.
        term_files : list of str
            All available terms files.
        gene_files : list of str
            All available gene files.
        """

        # Get lists of files from data directories
        osc_files = clean_file_list(os.listdir(os.path.join(self.maps_path, 'Oscs')), 'osc')
        slope_files = clean_file_list(os.listdir(os.path.join(self.maps_path, 'Slopes')), 'slope')
        gene_files = clean_file_list(os.listdir(os.path.join(self.maps_path, 'Genes')), 'gene')
        term_files = clean_file_list(os.listdir(os.path.join(self.maps_path, 'Terms')), 'terms')

        # If asked for, print out lists of files
        if verbose:
            print('Oscillation Files:\n', '\n'.join(osc_files), '\n')
            print('Slope Files:\n', '\n'.join(slope_files), '\n')
            print('Terms Files:\n', '\n'.join(term_files), '\n')
            print('Genes Files:\n', '\n'.join(gene_files), '\n')

        # If asked for, return lists of files
        if return_files:
            return osc_files, slope_files, term_files, gene_files

#################################################################################################
############################## OMEGAMAPPIN - CORE - DB - FUNCTIONS ##############################
#################################################################################################

def make_file_directory_internal(base_path):
    """Generates the database folder structure for internal data."""

    # Corrs Data
    os.mkdir(os.path.join(base_path, 'Corrs'))
    cor_data = ['Genes', 'Terms']
    cor_data_type = ['csv', 'npz']

    for cor_dat in cor_data:
        os.mkdir(os.path.join(base_path, 'Corrs', cor_dat))

        for cor_dat_type in cor_data_type:
            os.mkdir(os.path.join(base_path, 'Corrs', cor_dat, cor_dat_type))

    # Maps Data
    os.mkdir(os.path.join(base_path, 'Maps'))
    maps_data = ['Genes', 'Oscs', 'Slopes', 'Terms', 'Anat', 'Scouts']
    for maps_dat in maps_data:
        os.mkdir(os.path.join(base_path, 'Maps', maps_dat))

    # Processed Data
    os.mkdir(os.path.join(base_path, 'Processed'))
    proc_data = ['meg', 'maps']
    for proc_dat in proc_data:
        os.mkdir(os.path.join(base_path, 'Processed', proc_dat))


def make_file_directory_external(base_path):
    """Generates the database folder structure for external data."""

    # MEG Data
    os.mkdir(os.path.join(base_path, 'MEG'))

    meg_data = ['FOOF', 'PSDs', 'Viz']
    for meg_dat in meg_data:
        os.mkdir(os.path.join(base_path, 'MEG', meg_dat))


def make_test_file_directory_other(base_path):
    """Generates the other folder for a test directory."""

    # Other directory
    os.mkdir(os.path.join(base_path, 'csvs'))

#################################################################################
################## OMEGAMAPPIN - CORE - DB - PRIVATE FUNCTIONS ##################
#################################################################################

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
