"""Database structure object for the OM project."""

import os

from om.core.utils import clean_file_list, get_sub_nums

###################################################################################################
###################################################################################################

class OMDB(object):
    """Class to hold database information for MEG project.

    Attributes
    ----------
    internal_path : str
        Base path to all internal data.
    external_path : str
        Base path to all external data.
    maps_path : str
        Path to Maps data.
    corrs_path : str
        Path to Corrs data.
    meg_save_path : str
        Path to save md data.
    maps_save_path : str
        Path to save mc data.
    meg_path : str
        Path to MEG data.
    psd_path : str
        Path to PSD data.
    fooof_path : str
        Path to FOOOF data.
    viz_path : str
        Path to vizualization data.
    """

    def __init__(self, auto_gen=True):
        """Initialize OMDB object."""

        # Initialize base paths
        self.internal_path = ("/Users/tom/Documents/Research/"
                              "1-Projects/MEGmapping/2-Data/OMData/")
        self.external_path = ("/Users/tom/Documents/Research/"
                              "1-Projects/MEGmapping/2-Data/ExternalData/")

        # Initialize all internal paths
        self.maps_path = str()
        self.corrs_path = str()
        self.save_path = str()

        # Initialize all external paths
        self.meg_path = str()
        self.psd_path = str()
        self.fooof_path = str()
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
        self.save_path = os.path.join(self.internal_path, 'Processed')

        # Set up external data paths
        self.meg_path = os.path.join(self.external_path, 'MEG')
        self.psd_path = os.path.join(self.meg_path, 'PSDs')
        self.fooof_path = os.path.join(self.meg_path, 'FOOOF')
        self.viz_path = os.path.join(self.meg_path, 'Viz')


    def check_dat_files(self, dat_type, dat_source='both', save_type='', verbose=False):
        """Checks what data files are available.

        Parameters
        ----------
        dat_type : {'PSD', 'fooof'}
            Which data type to check files for.
        dat_source : {'OMEGA', 'HCP', 'both'}
            Which database to check files for.
        verbose : boolean, optional (default = False)
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
            f_l = 'last'
        elif dat_type is 'fooof':
            dat_path = self.fooof_path
            word = 'fooof'
            f_l = 'first'
        else:
            raise ValueError('Data type not understood.')

        # If looking for a particular database, find file, get subject numbers and source
        if dat_source is not 'both':
            sub_nums = _check_files(os.path.join(dat_path, dat_source), word, f_l)
            source = [dat_source] * len(sub_nums)

        # If looking across both databases, get info from each database and then combine
        else:
            sub_nums_omega = _check_files(os.path.join(dat_path, 'OMEGA'), word, f_l)
            n_omega = len(sub_nums_omega)

            sub_nums_hcp = _check_files(os.path.join(dat_path, 'HCP'), word, f_l)
            n_hcp = len(sub_nums_hcp)

            sub_nums = sub_nums_omega + sub_nums_hcp
            source = (['OMEGA'] * n_omega) + (['HCP'] * n_hcp)

        # If requested, print out the list of subject numbers
        if verbose:
            print('\nNumber of Subjects available: ' + str(len(sub_nums)) + '\n')
            print('Subject numbers with FOOOF data available: \n' + str(sub_nums) + '\n')

        return sub_nums, source


    def check_res_files(self, dat_type, verbose=True):
        """Checks what result files are available.

        Parameters
        ----------
        res_type : {'meg', 'maps'}
            Which data type to check files for.
        verbose : boolean, optional (default = True)
            Whether to print out information during run.

        Returns
        -------
        files : list of str
            A list of all the available files.
        """

        # Get files
        files = os.listdir(os.path.join(self.save_path, dat_type))
        files = clean_file_list(files, dat_type)

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

def check_db(db):
    """Checks if OMDB object is initialized. If not, returns OMDB() object.

    Parameters
    ----------
    db : {OMDB() object, None}
        Can be None, or a OMDB() object.

    Returns
    -------
    db : OMDB() object
        Database object for OM project.
    """

    if not db:
        db = OMDB()

    return db

#########################################################################################
########################## CREATE DATABASE STRUCTURE FUNCTIONS ##########################
#########################################################################################

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

    meg_data = ['FOOOF', 'PSDs', 'Viz']
    for meg_dat in meg_data:
        os.mkdir(os.path.join(base_path, 'MEG', meg_dat))


def make_test_file_directory_other(base_path):
    """Generates the other folder for a test directory."""

    # Other directory
    os.mkdir(os.path.join(base_path, 'csvs'))

###################################################################################################
###################################################################################################

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
