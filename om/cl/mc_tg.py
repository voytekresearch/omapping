"""MODULE DOCSTRING - TO FILL IN

"""

from __future__ import print_function, division
import os
import csv
import numpy as np
import pandas as pd
import scipy.stats.stats as sps
from ipyparallel import Client
from ipyparallel.util import interactive

# Import custom om code
from om.gen import OMDB, DataNotComputedError, UnknownDataTypeError, clean_file_list
from om.cl.mc_base import MapCompBase

######################################################################################
########################## OMEGAMAPPIN - CL_MC_TG - CLASSES ##########################
######################################################################################

class MapCompTG(MapCompBase):
    """DOCSTRING"""

    def __init__(self, db):
        """

        Parameters
        ----------
        db : xx
            xx
        """

        # Inherit from MapCompBase() class
        MapCompBase.__init__(self, db)

        # Import the vectors of gene & term names
        self.term_names = _get_map_names('00-ns_terms.csv', self.db.maps_terms_path)
        self.gene_names = _get_map_names('00-real_gene_names.csv', self.db.maps_genes_path)

        # Get number of terms and genes used
        self.n_terms = len(self.term_names)
        self.n_genes = len(self.gene_names)

        # Initialize variable to store the term maps
        self.term_maps = np.array([])

        # Initialize variable to store the gene maps
        self.gene_maps = np.array([])
        self.gene_subj = str()

        # Initialize a dictionary to store all the R-value results from spatial correlations
        self.corrs = dict({'Terms': dict(), 'Genes': dict()})

        # Initialize a dictionary to store all the p-value results from spatial correlations
        self.p_vals = dict({'Terms': dict(), 'Genes': dict()})

        # Initialize booleans that keep track of what is loaded
        self.terms_loaded = False
        self.genes_loaded = False


    def load_gene_maps(self, subject):
        """Load the spatial maps of gene data.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        subject : str
            Which subject of gene data to load. Of the form 'sub#'
        """

        # Check if gene data already loaded - if so, unload
        if self.genes_loaded:
            print('Unloading previously loaded genes.')
            self.unload_data('Genes')

        # Set subject marker
        self.gene_subj = subject

        # Get the list of gene file names for given subject
        genes_file_names = _get_gene_files(subject)

        # Check how many files there are
        n_files = len(genes_file_names)

        # Loop through all files given
        for ind, genes_csv in enumerate(genes_file_names):

            # Print loading status
            print('Loading file #', str(ind+1), ' of ', str(n_files))

            # If first file, initialize as first part of the gene map
            if ind == 0:
                self.gene_maps = pd.read_csv(genes_csv, header=None)

            # For all subsequent files, concatenate to end of gene map
            else:
                temp_df = pd.read_csv(genes_csv, header=None)
                self.gene_maps = pd.concat([self.gene_maps, temp_df])

        # Print status
        print('All files loaded!')

        # Check that term data loaded matches number of term names
        [n_vert, n_genes] = self.gene_maps.shape
        if n_genes != self.n_genes:
            print('NUMBER OF GENES DOES NOT MATCH')

        # Update boolean that genes are loaded
        self.genes_loaded = True


    def load_term_maps(self, term_file_name):
        """Load the spatial maps of term data.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        term_file_name : str
            File name of term data file.
        """

        # Get full path for the csv file
        terms_csv = os.path.join(self.db.maps_terms_path, term_file_name)

        # Load the terms map
        self.term_maps = pd.read_csv(terms_csv, header=None)

        # Check that term data loaded matches number of term names
        [n_vert, n_terms] = self.term_maps.shape
        if n_terms != self.n_terms:
            print('NUMBER OF TERMS DOES NOT MATCH')

        # Update boolean that terms are loaded
        self.terms_loaded = True


    def calc_corrs(self, dat_type, meg_dat, method='linear'):
        """Calculate correlations between spatial maps.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        dat_type : str
            Type of data to correlate with meg data.
                Options: {'Terms', 'Genes'}
        meg_dat : str
            Specific type of meg data to correlate.
                osc_band or 'Slopes' only
        method : {'linear', 'parallel'}
            Run method (linear or parallel) to use.
        """

        # Check with data type and set data accordingly
        if dat_type is 'Terms':
            n_comps = self.n_terms
            dat_df = self.term_maps
        elif dat_type is 'Genes':
            n_comps = self.n_genes
            dat_df = self.gene_maps
        else:
            raise UnknownDataTypeError('Data Type not understood.')

        # Get the specified meg map
        if meg_dat is 'Slopes':

            # Check that slopes are loaded
            if not self.slopes_loaded:
                raise DataNotComputedError('Slope data has not been loaded.')
            meg_map = self.slope_map[meg_dat]

        else:

            # Check that oscillation data is loaded
            if not self.oscs_loaded:
                raise DataNotComputedError('Oscillation data has not been loaded.')
            try:
                meg_map = self.meg_maps[meg_dat]
            except KeyError:
                raise UnknownDataTypeError('MEG Data not understood.')

        # Initialize dictionaries to store correlation data
        if not self.corrs[dat_type]:
            self.corrs[dat_type] = _init_stat_dict(self.bands)
            self.p_vals[dat_type] = _init_stat_dict(self.bands)

        # Print out status
        print('Calculating corrs between', str(dat_type), 'and', str(meg_dat))

        # Run linearly
        if method is 'linear':

            # Initialize vectors to store correlation results
            corr_vals = np.zeros([n_comps, 1])
            p_vals = np.zeros([n_comps, 1])

            # Loop through all comparisons to run
            for comp in range(0, n_comps):

                # Pull out specific data (single term or gene)
                dat = np.array(dat_df.ix[:, comp])

                # Get inds of data that contains numbers
                inds_non_nan = [i for i in range(len(dat)) if not np.isnan(dat[i])]

                # Calculate correlation between data and meg map
                [corr_vals[comp], p_vals[comp]] = sps.pearsonr(
                    dat[inds_non_nan], meg_map[inds_non_nan])

        # Run in parallel
        elif method is 'parallel':

            # Initialize client & gather workers
            client = Client()
            view = client[:]

            # Import required libraries for each worker
            with view.sync_imports():
                import numpy
                from scipy.stats.stats import pearsonr

            # Send data to workers
            view['meg_map'] = meg_map

            # Turn data into a list
            dat_list = _make_list(dat_df)

            # Map and get results
            corr_map = view.map(_run_corr, dat_list)
            results = corr_map.get()

            # Pull out results
            [corr_vals, p_vals] = _pull_out_results(results)

        # Save correlations results to MapComp object
        self.corrs[dat_type][meg_dat] = corr_vals
        self.p_vals[dat_type][meg_dat] = p_vals


    def check_corrs(self, dat_type, meg_dat, n_check=20, top=True):
        """Check (print out) highest or lowest oscillations.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        dat_type : {'Terms', 'Genes'}
            Data type (Terms or Genes) of corrs to check.
        meg_dat : str
            Specific MEG data of corrs to check.
        n_check : int, optional (default = 20)
            Number of correlations to check.
        top : boolean, optional (default = True)
            Get Top (True) or Bottom (False) set of correlations.
        """

        # Check which type of data and set names accordingly
        if dat_type is 'Terms':
            names = self.term_names
            la_str = ' <30'
        elif dat_type is 'Genes':
            names = self.gene_names
            la_str = ' <55'
        else:
            raise UnknownDataTypeError('Data type not understood.')

        # Check that asked for correlations have been computed
        if not self.corrs[dat_type]:
            raise DataNotComputedError('No correlations calculated for requested data type.')
        elif len(self.corrs[dat_type][meg_dat]) == 0:
            raise DataNotComputedError('Requested meg data correlation not calculated.')

        # Get R and p values of specified correlations
        meg_corr = np.squeeze(self.corrs[dat_type][meg_dat])
        meg_p = np.squeeze(self.p_vals[dat_type][meg_dat])

        # Sort the corr vector
        inds_max = np.argsort(meg_corr, axis=0)

        # Print Header Rows
        print("\n\nCorrelations for ", str(dat_type), " &  ", str(meg_dat), ': \n')
        print('# \t', format(dat_type, la_str), '\t R-Vals \t P-vals \n')

        # Print out the top group (highest positive correlations)
        if top:
            for i in range(1, n_check+1):
                ind = int(inds_max[-i])
                print(str(i), '\t', format(names[ind][:50], la_str), '\t',
                      format(meg_corr[ind], '1.5f'), '\t', format(meg_p[ind], '1.4e'))

        # Print out the bottom group (highest negative correlations)
        else:
            for i in range(0, n_check):
                ind = int(inds_max[i])
                print(str(i), '\t', format(names[ind][:50], la_str), '\t',
                      format(meg_corr[ind], '1.5f'), '\t', format(meg_p[ind], '1.4e'))


    def unload_data(self, dat_type):
        """Unload specified data from MapComp object.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        dat_type : {'Terms', 'Genes'}
            Indicator of which set of data to unload from object.
        """

        # Unload Term data
        if dat_type is 'Terms':

            # Check if terms are currently loaded. Return if not.
            if not self.terms_loaded:
                raise DataNotComputedError('Terms not loaded - can not unload.')

            # Unload terms by resetting map variable as empty
            self.term_maps = np.array([])

            # Update boolean that terms are not loaded
            self.terms_loaded = False

        # Unload Gene data
        elif dat_type is 'Genes':

            # Check if genes are currently loaded. Return if not.
            if not self.genes_loaded:
                raise DataNotComputedError('Terms not loaded - can not unload.')

            # Unload genes by resetting map variable as empty
            self.gene_maps = np.array([])
            self.gene_subj = str()

            # Update boolean that genes are not loaded
            self.genes_loaded = False

        # Otherwise, data type was not understood
        else:
            raise UnknownDataTypeError('Data type not understood.')


    def save_corrs(self, dat_type, meg_dat, save_title, save_as_npz=True, save_as_csv=True):
        """Save out the correlation results.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        dat_type : {'Terms', 'Genes'}
            Data type to save corrs for.
        meg_dat : str
            MEG data to save corrs for.
        save_title : str
            String to attach to the front of the save file names.
        save_as_npz : boolean, optional (default = True)
            Whether to save an npz file.
        save_as_csv : boolean, optional (default = True)
            Whether to save a csv file.
        """

        # Check which type of data and set names, filenames & save paths accordingly
        if dat_type is 'Terms':
            names = self.term_names
            file_name = save_title + '_Corrs_' + dat_type + '_' + meg_dat
            save_path = os.path.join(self.db.corrs_path, dat_type)
            sub_name = ''
        elif dat_type is 'Genes':
            names = self.gene_names
            file_name = save_title + '_' + self.gene_subj + '_Corrs_' + dat_type + '_' + meg_dat
            save_path = os.path.join(self.db.corrs_path, dat_type)
            sub_name = self.gene_subj
        else:
            raise UnknownDataTypeError('Data type not understood.')

        # Check that asked for correlations have been computed
        if not self.corrs[dat_type]:
            raise DataNotComputedError('No correlations calculated for requested data type.')
        elif len(self.corrs[dat_type][meg_dat]) == 0:
            raise DataNotComputedError('Requested meg data correlation not calculated.')

        # Get the correlation data of interest
        meg_corrs = np.squeeze(self.corrs[dat_type][meg_dat])
        meg_p_vals = np.squeeze(self.p_vals[dat_type][meg_dat])

        # Save a numpy npz file
        if save_as_npz:

            # Set up specific name/path for npz file
            npz_path = os.path.join(save_path, 'npz', sub_name, '')
            outfile_npz = npz_path + file_name + '.npz'

            # Save out npz file
            np.savez(outfile_npz, dat_type=dat_type, meg_dat=meg_dat, names=names,
                     meg_corrs=meg_corrs, meg_p_vals=meg_p_vals)

        # Save a csv file
        if save_as_csv:

            # Set up specific name/path for csv file
            csv_path = os.path.join(save_path, 'csv', sub_name, '')
            outfile_csv = csv_path + file_name + '.csv'

            # Open the csv file
            csv_file = open(outfile_csv, 'w')

            # Write Header
            csv_file.write('Name, R-Value, P-Value' + '\n')

            # Get number of rows of data to save out
            n_rows = len(names)

            # Sort data to save out in sorted order
            inds_max = np.argsort(meg_corrs, axis=0)

            # Loop through each row of data, saving out to file
            for i in range(1, n_rows+1):

                # Get index of next data to save
                ind = int(inds_max[-i])

                # Set data as string to print out to csv file
                out_list = [names[ind].replace(',', ' '), str(meg_corrs[ind]), str(meg_p_vals[ind])]
                out_str = ",".join(out_list)
                csv_file.write(out_str + '\n')

            # Close the csv file
            csv_file.close()

###############################################################################################
######################### OMEGAMAPPIN - CL_MC_TG - FUNCTIONS (PUBLIC) #########################
###############################################################################################

def calc_avg_gene_map(subj_list, file_title):
    """

    Parameters
    ----------
    subj_list : list of str
        List of subject numbers to average together
    file_title : xx
        xx

    Outputs
    -------
    xx : xx
        xx
    """

    # Get OMDB object
    db = OMDB()

    # Check how many subjects to average over
    n_subjs = len(subj_list)

    #
    in_files_path = []
    for subj in subj_list:
        in_files_path.append(_get_gene_files(subj))

    # Loop through three file parts
    for part in range(3):

        # Set up output file
        out_file = file_title + '_genes_average_r10_' + str(part+1) + 'of3.csv'
        out_file_path = os.path.join(db.maps_genes_path, 'avg_gene_estimations', out_file)

        # Get current set of input files
        cur_part_in_files = []
        for subj in range(n_subjs):
            cur_part_in_files.append(in_files_path[subj][part])

        # Save out an average csv file from input files
        _avg_csv_files(cur_part_in_files, out_file_path)


################################################################################################
######################### OMEGAMAPPIN - CL_MC_TG - FUNCTIONS (PRIVATE) #########################
################################################################################################

def _get_map_names(names_file, path):
    """Get the map names from a given file.

    Parameters
    ----------
    names_file : str
        File name to pull the map names from.
    path : str
        Path of where the file is.

    Returns
    -------
    names : list of str
        A list of the map names
    """

    # Get path to csv file
    csv_path = os.path.join(path, names_file)

    # Open csv file
    with open(csv_path, 'rb') as f_name:
        reader = csv.reader(f_name, delimiter=',')

        # Get list of names from first row in csv
        names = list(reader)[0]

    return names


def _get_gene_files(subj):
    """Returns full file paths for all gene files for a given subject.

    Parameters
    ----------
    subj : str
        Which subject to get files for

    Outputs
    -------
    file_names_path : list of str
        A list of full file names for all gene files for given subject
    """

    # Get OMDB object
    db = OMDB()
    #maps_genes_path = os.path.join(db.maps_path, 'Genes')

    # Make var for the name of the folder of gene data
    subj_folder = subj + '_gene_estimations'

    # Get a list of all the files in the folder
    file_names = clean_file_list(os.listdir(os.path.join(
        db.maps_genes_path, subj_folder)), 'r10')

    # Check how many files there are
    n_files = len(file_names)

    # Make a list of the full file names, including full path
    file_names_path = []
    for cur_file in range(n_files):
        file_names_path.append(os.path.join(db.maps_genes_path, subj_folder, file_names[cur_file]))

    return file_names_path


def _avg_csv_files(f_in, f_out, avg='mean'):
    """

    Note: This function assumes csv files with a constant number of rows
        and columns, the same for all input files. Will fail, perhaps
        silently if this is not the case.

    Parameters
    ----------
    f_in : ?
        xx
    f_out : ?
        xx
    avg : {'mean', 'median'}
        xx
    """

    # Open out file object
    out_file = open(f_out, 'wb')
    out_writer = csv.writer(out_file)

    # Check how many input files there are
    n_in = len(f_in)

    # Create input file objects
    in_files = []
    in_readers = []
    for i in range(n_in):
        in_files.append(open(f_in[i]))
        in_readers.append(csv.reader(in_files[i]))

    # Initialize number of columns as false
    n_col = False

    # Loop through each line of
    for row in in_readers[0]:

        # If unknown, check the number of columns
        if not n_col:
            n_col = len(row)

        # Initialize a temporary array to store
        temp = np.zeros([n_in, n_col])

        # Add first row to
        temp[0, :] = np.array([float(i) for i in row])

        #
        for f_ind in range(1, n_in):

            # Load row of data into the temporary array
            temp[f_ind, :] = np.array([float(i) for i in in_readers[f_ind].next()])

        # Take average
        if avg is 'mean':
            avg_dat = np.nanmean(temp, 0)
        elif avg is 'median':
            avg_dat = np.nanmedian(temp, 0)

        # Write out line to average csv file
        out_writer.writerow(avg_dat.tolist())

    # Close out all files
    for i in range(n_in):
        in_files[i].close()
    out_file.close()


def _init_stat_dict(bands):
    """Initialize a dictionary to store stat data for inter-data correlations.

    Parameters
    ----------
    bands : list of str
        Oscillation bands to initialize.
    """

    # Initialize dictionary to return
    out = dict()

    # Add to dictionary for each band
    for band in bands:
        out[band] = np.array([])

    # Add a field for slope correlations
    out['Slopes'] = np.array([])

    return out


def _make_list(dat_df):
    """Turns a dataframe into a list.

    Parameters
    ----------
    dat_df : DataFrame
        Input data, to be converted into a list.
    """

    # Get size of the data to
    [n_verts, n_dat] = dat_df.shape

    # Initialize list to return
    out_list = list()

    # Pull each data column into a list entry
    for i in range(n_dat):
        out_list.append(np.array(dat_df.ix[:, i]))

    return out_list


def _pull_out_results(dat_in):
    """Given a list of correlation results, pulls them out into arrays.

    Parameters
    ----------
    dat_in : list of tuple
        A list of correlation results. Each tuple is (R-val, p-val).
    """

    # Check length of data
    n_dat = len(dat_in)

    # Initializ vectors
    out_1 = np.zeros([n_dat, 1])
    out_2 = np.zeros([n_dat, 1])

    # Loop through and pull out data
    for i in range(n_dat):
        out_1[i] = dat_in[i][0]
        out_2[i] = dat_in[i][1]

    return out_1, out_2


@interactive
def _run_corr(dat):
    """Run correlation between maps. Used for parallel runs.

    Parameters
    ----------
    dat : 1d array
        An array of map data to be compared to projected meg map.

    Note:
    - meg_map has to be projected to workers.
    - numpy and pearsonr have to be imported on workers.
    """

    # Get inds of data that contains numbers
    inds_non_nan = [i for i in range(len(dat)) if not numpy.isnan(dat[i])]

    # Calculate corr between data and MEG map
    [corr_vals, p_vals] = pearsonr(dat[inds_non_nan], meg_map[inds_non_nan])

    return (corr_vals, p_vals)
