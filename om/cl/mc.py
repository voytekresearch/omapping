from __future__ import print_function, division
import os
import csv
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.stats.stats as sps
import matplotlib.pyplot as plt
from ipyparallel import Client
from ipyparallel.util import interactive

# Import custom om code
from om.gen import *

# TODO: Update anything in here relevant to the update to Osc_Dict().
#   Anything with hard coded oscillation bands need fixing. 

###########################################################################################
############################ OMEGAMAPPIN - MAP COMPARE CLASSES ############################
###########################################################################################

class MapComp():
    """Class for storing and comparing spatial topographies."""

    def __init__(self, db):

        # Pull out needed paths from OMDB object
        self.project_path = db.project_path
        self.maps_path = db.maps_path
        self.corrs_path = db.corrs_path

        # Set specific paths for different data types
        self.oscs_path = os.path.join(self.maps_path, 'Oscs')
        self.slopes_path = os.path.join(self.maps_path, 'Slopes')
        self.terms_path = os.path.join(self.maps_path, 'Terms')
        self.genes_path = os.path.join(self.maps_path, 'Genes')

        # Import the vectors of gene & term names
        self.term_names = _get_map_names('ns_terms.csv', self.terms_path)
        self.gene_names = _get_map_names('real_gene_names.csv', self.genes_path)

        # Get number of terms and genes used
        self.n_terms = len(self.term_names)
        self.n_genes = len(self.gene_names)

        # Initialize a dictionary to store maps of meg data (oscillations & slopes)
        self.meg_maps = dict()
        #self.meg_maps = _init_meg_map_dict()

        # Initialize variable to store the term maps
        self.term_maps = np.array([])

        # Initialize variable to store the gene maps
        self.gene_maps = np.array([])
        self.gene_subj = str()

        # Initialize a dictionary to store all the R-value results from spatial correlations
        self.corrs = _init_stat_dict()

        # Initialize a dictionary to store all the p-value results from spatial correlations
        self.p_vals = _init_stat_dict()

        # Initialize booleans that keep track of what is loaded
        self.oscs_loaded = False
        self.slopes_loaded = False
        self.terms_loaded = False
        self.genes_loaded = False


    def check_files(self, print_files=True, return_files=False):
        """Gets the list of files in the map directories. Can return and/or print.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        print_files : boolean, optional
            Whether or not to print out available file names.
        return_files : boolean, optional
            Whether or not to return lists of filenames.
        """

        # Get lists of files from data directories
        osc_files = clean_file_list(os.listdir(self.oscs_path), 'osc')
        slope_files = clean_file_list(os.listdir(self.slopes_path), 'slope')
        gene_files = clean_file_list(os.listdir(self.genes_path), 'gene')
        term_files = clean_file_list(os.listdir(self.terms_path), 'terms')

        # If asked for, print out lists of files
        if print_files:
            print('Oscillation Files:\n', '\n'.join(osc_files), '\n')
            print('Slope Files:\n', '\n'.join(slope_files), '\n')
            print('Terms Files:\n', '\n'.join(term_files), '\n')
            print('Genes Files:\n', '\n'.join(gene_files), '\n')

        # If asked for, return lists of files
        if return_files:
            return osc_files, slope_files, term_files, gene_files


    def load_meg_maps(self, osc_file=None, slope_file=None):
        """Load the spatial maps of MEG data (oscillations and slopes).

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        osc_file : str, optional
            File path to the file with oscillations data.
        slope_file : str, optional
            File path to the file with slope data.
        """

        # Initialize the var to store meg map data
        if slope_file is not None:
            self.meg_maps = _init_meg_map_dict()
        else:
            self.meg_maps = _init_meg_map_dict(slopes=False)

        # If a filename is provided, load oscillation data
        if osc_file is not None:

            # Get the full path for the file name
            oscs_map_file = os.path.join(self.oscs_path, osc_file + '.npz')

            # Load osc file data
            with np.load(oscs_map_file) as data:

                # Load osc maps into the meg_map dictionary
                self.meg_maps['Theta']      = data['osc_score_theta']
                self.meg_maps['Alpha']      = data['osc_score_alpha']
                self.meg_maps['Beta']       = data['osc_score_beta']
                self.meg_maps['LowGamma']   = data['osc_score_lowgamma']

            # Update boolean that oscs are loaded
            self.oscs_loaded = True

        # If a filename is provided, load slope data
        if slope_file is not None:

            # Get the full path for the file name
            slopes_map_file = os.path.join(self.slopes_path, slope_file + '.npz')

            # Load slope file data
            with np.load(slopes_map_file) as data:

                # Load the slope map into the meg_map dictionary
                self.meg_maps['Slopes']     = data['chis']

            # Update boolean that slopes are loaded
            self.slopes_loaded = True


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
            self.unload_data('Genes')

        # Set subject marker
        self.gene_subj = subject

        # Make string for folder name of subject gene data directory
        subj_str = subject + '_gene_estimations'

        # Get list of files available for requested subject
        genes_file_names = clean_file_list(os.listdir(os.path.join(self.genes_path, subj_str)), 'r10')

        # If one file given, load this as the gene map
        if len(genes_file_names) == 1:
            genes_csv = os.path.join(self.genes_path, subj_str, genes_file_names[0])
            self.gene_maps = pd.read_csv(genes_file_names[0], header=None)

        # If multiple files, load them all and concatenate
        else:

            # Check how many files there are
            nFiles = len(genes_file_names)

            # Loop through all files given
            for i in range(0, nFiles):

                # Print loading status
                print('Loading file #', str(i+1), ' of ', str(nFiles))

                # If first file, initialize as first part of the gene map
                if i == 0:
                    genes_csv = os.path.join(self.genes_path, subj_str, genes_file_names[0])
                    self.gene_maps = pd.read_csv(genes_csv, header=None)

                # For all subsequent files, concatenate to end of gene map
                else:
                    genes_csv = os.path.join(self.genes_path, subj_str, genes_file_names[i])
                    temp_df = pd.read_csv(genes_csv, header=None)
                    self.gene_maps = pd.concat([self.gene_maps, temp_df])

        # Print status
        print('All files loaded!')

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
        terms_csv = os.path.join(self.terms_path, term_file_name)

        # Load the terms map
        self.term_maps = pd.read_csv(terms_csv, header=None)

        # Update boolean that terms are loaded
        self.terms_loaded = True


    def calc_corrs(self, dat_type, meg_dat, method='linear'):
        """Calculate correlations between spatial maps.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        dat_type : str
            Type of data to correlate with meg data
                'Terms' or 'Genes' only
        meg_dat : str
            Specific type of meg data to correlate
                osc_band or 'Slopes' only
        method : str
            Run method (linear or parallel) to use.
                Options: {'linear', 'parallel'}. Default: 'linear'.
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
        try:
            meg_map = self.meg_maps[meg_dat]
        except KeyError:
            raise UnknownDataTypeError('MEG Data not understood.')

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
                [corr_vals[comp], p_vals[comp]] = sps.pearsonr(dat[inds_non_nan], meg_map[inds_non_nan])

        # Run in parallel
        elif method is 'parallel':

            # Initialize client & gather workers
            c = Client()
            view = c[:]

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
        self.corrs[dat_type + meg_dat] = corr_vals
        self.p_vals[dat_type + meg_dat] = p_vals


    def check_corrs(self, dat_type, meg_dat, n_check=20, top=True):
        """Check (print out) highest or lowest oscillations.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        dat_type : str
            Data type (Terms or Genes) of corrs to check.
        meg_dat : str
            Specific MEG data of corrs to check.
        n_check : int, optional
            Number of correlations to check.
        top : boolean, optional.
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
        #if not self.corrs[dat_type + meg_dat]:
        #    print("Those correlations not calculated. Quitting.")
        #    return

        # Get R and p values of specified correlations
        meg_corr = np.squeeze(self.corrs[dat_type + meg_dat])
        meg_p = np.squeeze(self.p_vals[dat_type + meg_dat])

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
        dat_type : str
            Indicator of which set of data to unload from object.
                Options: 'Terms', 'Genes'
        """

        # Unload Term data
        if dat_type is 'Terms':

            # Check if terms are currently loaded. Return if not.
            if not self.terms_loaded:
                print("Terms are not loaded - can not unload.")

            # Unload terms by resetting map variable as empty
            self.term_maps = np.array([])

            # Update boolean that terms are not loaded
            self.terms_loaded = False

        # Unload Gene data
        elif dat_type is 'Genes':

            # Check if genes are currently loaded. Return if not.
            if not self.genes_loaded:
                print("Genes are not loaded - can not unload.")

            # Unload genes by resetting map variable as empty
            self.gene_maps = np.array([])
            self.gene_subj = str()

            # Update boolean that genes are not loaded
            self.genes_loaded = False

        # Otherwise, data type was not understood
        else:
            raise UnknownDataTypeError('Data type not understood.')


    def save_corrs(self, dat_type, meg_dat, save_as_npz=True, save_as_csv=True):
        """Save out the correlation results.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        dat_type : str
            Data type to save corrs for.
                Options: {'Terms', 'Genes'}
        meg_dat : str
            MEG data to save corrs for.
                Options: {'Theta', 'Alpha', 'Beta', 'LowGamma', 'Slopes'}
        save_as_npz : boolean, optional
            Whether to save an npz file. Default is True.
        save_as_csv : boolean, optional
            Whether to save a csv file. Default is True.
        """

        # Check which type of data and set names, filenames & save paths accordingly
        if dat_type is 'Terms':
            names = self.term_names
            file_name = 'Corrs_' + dat_type + '_' + meg_dat
            save_path = os.path.join(self.corrs_path, dat_type)
            sub_name = ''
        elif dat_type is 'Genes':
            names = self.gene_names
            file_name = self.gene_subj + '_Corrs_' + dat_type + '_' + meg_dat
            save_path = os.path.join(self.corrs_path, dat_type)
            sub_name = self.gene_subj
        else:
            raise UnknownDataTypeError('Data type not understood.')

        # Check that asked for correlations have been computed
        #if not self.corrs[dat_type + meg_dat]:
        #    print("Those correlations not calculated. Quitting.")
        #    return

        # Get the correlation data of interest
        meg_corrs = np.squeeze(self.corrs[dat_type + meg_dat])
        meg_p_vals = np.squeeze(self.p_vals[dat_type + meg_dat])

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
                out_list = ",".join(out_list)
                csv_file.write(out_list + '\n')

            # Close the csv file
            csv_file.close()


class MapCompROI(MapComp):
    """Class for storing and comparing spatial topographies in ROIs."""

    def __init__(self, db):

        # Inherit from MapComp() class
        MapComp.__init__(self, db)

        # Initialize var to store number of ROIs
        self.nROIs = int()

        # Add path for anatomy data
        self.anat_path = os.path.join(self.maps_path, 'Anat')

        # Add vars to save ROI data from anat data
        self.anat_roi_names = list()
        self.anat_roi_lrs = list()
        self.anat_nROIs = int
        self.anat_con = np.ndarray(0)
        self.anat_type = str()

        # Add vars to save ROI data from MEG data
        self.elec_roi_names = list()
        self.elec_roi_lrs = list()
        self.elec_nROIs = int
        self.elec_roi_verts = list()

        # Initialize list to store ROI labels
        self.roi_labels = list()
        self.roi_verts = list()
        self.roi_lr = list()
        self.nROIs = int

        # Initialize var to store MEG ROI data
        self.meg_ROI_maps = dict()

        # Initialize var to store MEG connectivity data
        self.meg_con = dict()

        # Initialize var to store meg stats data
        self.meg_stats = dict()

        # Add boolean for whether anat data is loaded
        self.anat_loaded = False
        self.elec_loaded = False
        self.rois_aligned = False
        self.meg_ROIs = False


    def load_anat_maps(self, anat_file_name, anat_type):
        """Load the spatial maps of anatomilcal data.

        Parameters
        ----------
        self : MapCompROI() object
            Object for storing and comparing map data, in ROI format.
        anat_file_name : str
            File name of anat data file.
        anat_type : str
            Indicates which type of anat data is loaded.
        """

        # Get full path for the anat mat file
        anat_mat_file = os.path.join(self.anat_path, anat_file_name)

        # Load the anat data
        dat = sio.loadmat(anat_mat_file)

        # Pull out data from mat dictionary
        roi_names = dat['roi_labels'].tolist()
        self.anat_con = dat['connectivity']

        # Get number of ROIs
        self.anat_nROIs = len(roi_names)

        # Loop through and fix roi labels
        for r in range(0, self.anat_nROIs):
            self.anat_roi_names.append(str(roi_names[r][0][0]))

        # Update which type of anatomy data is loaded
        self.anat_type = anat_type

        # Update boolean that anat data is loaded
        self.anat_loaded = True


    def load_elec_rois(self, roi_file_name):
        """Load the ROI file for the MEG data.

        Parameters
        ----------
        self : MapCompROI() object
            Object for storing and comparing map data, in ROI format.
        roi_file_name : str
            File name (including path) to the matlab scout definition file.
        """

        # Load ROI dat from mat file
        dat = sio.loadmat(roi_file_name, struct_as_record=True)

        # Pull out data from mat file
        scouts = np.squeeze(dat['Scouts'])

        # Check how many ROIs there are
        self.elec_nROIs = len(scouts)

        # Initiliaze vars to store scout names and vertices
        sc_names = list()
        sc_verts = list()

        # Loop through, pull out names and verts into lists
        for i in range(0, self.elec_nROIs):
            sc = scouts[i]
            sc_verts.append(sc[0])
            sc_names.append(str(sc[3]))

            # Drop brackets in scout name
            sc_names[i] = sc_names[i][3:-2]

        # Attach data to object
        self.elec_roi_names = sc_names
        self.elec_roi_verts = sc_verts

        # Update boolean that elec data is loaded
        self.elec_loaded = True


    def align_rois(self):
        """Align ROIs used and names between anat and meg ROIs.

        NOTES
        -----
        xxxxxxxx
        """

        # Check if ROIs loaded - return if not
        if (not self.elec_roi_names) or (not self.anat_roi_names):
            print('One or Both ROIs not loaded! Cant proceed!')
            return

        # Elec L/Rs - standardize names
        for r in range(0, self.elec_nROIs):

            # Check if left side scout
            if self.elec_roi_names[r][-1] is 'L':
                self.elec_roi_lrs.append('L')

            # Check if right side scout
            elif self.elec_roi_names[r][-1] is 'R':
                self.elec_roi_lrs.append('R')

            else:
                pass

            # Drop the L/R from the names
            self.elec_roi_names[r] = self.elec_roi_names[r][:-2]

        # Anat L/Rs - standardize names
        for r in range(0, self.anat_nROIs):

            # Check if left side scout
            if self.anat_roi_names[r][0] is 'l':
                self.anat_roi_lrs.append('L')
                self.anat_roi_names[r] = self.anat_roi_names[r][5:]

            # Check if right side scout
            elif self.anat_roi_names[r][0] is 'r':
                self.anat_roi_lrs.append('R')
                self.anat_roi_names[r] = self.anat_roi_names[r][6:]

            else:
                pass

            # Drop the underscores
            self.anat_roi_names[r] = self.anat_roi_names[r].replace("_", "")

        # Loop through and line up scout names
        for i in range(0, self.anat_nROIs):

            # Grab current ROI from anat ROI list
            cur_roi = self.anat_roi_names[i]
            cur_lr = self.anat_roi_lrs[i]

            # Loop through elec ROIs to match up
            for j in range(0, self.elec_nROIs):

                # Check if current elec ROI matches current roi
                if self.elec_roi_names[j] == cur_roi:

                    # Check if L/R is right
                    if self.elec_roi_lrs[j] == cur_lr:

                        # Same side - add to overall list
                        self.roi_labels.append(cur_roi)
                        self.roi_lr.append(cur_lr)

                        # Add vertices to overall ROI list
                        self.roi_verts.append(self.elec_roi_verts[j])

        # Check how many ROIs there are
        self.nROIs = len(self.roi_labels)

        # Set boolean that ROIs have been aligned
        self.rois_aligned = True


    def conv_meg_rois(self):
        """Convert MEG data to ROIs.

        NOTES
        -----
        This XXXXX
        """

        # Initialize dict for current ROI data
        roi_meg_dat = _init_meg_map_dict(self.nROIs)

        # Loop through all ROIs
        for r in range(0, self.nROIs):

            # Add current ROI data to dict
            # Loop through all oscs
            for key in self.meg_maps.keys():

                #
                cur_verts = np.squeeze(self.roi_verts[r] - 1)
                nVerts = len(cur_verts)

                #
                temp_dat = self.meg_maps[key][cur_verts]

                #
                roi_meg_dat[key][r] = (sum(temp_dat) / nVerts)

        # Add the current ROI data to object
        self.meg_ROI_maps = roi_meg_dat

        # Update boolean that meg data has been converted to ROIs
        self.meg_ROIs = True


    def comp_meg_anat(self, section='all', print_out=True):
        """Compare anatomical connectivity to oscillation data.

        Parameters
        ----------
        self : MapCompROI() object
            Object for storing and comparing map data, in ROI format.
        section : str
            Which section of data to compare.
                Options: 'all', 'left', 'right'
        print_out : boolean, optional
            Whether to print out the stats results.
                Defaults to True
        """

        # Initialize the dictionary to store MEG connectivity data
        self.meg_con = _init_meg_map_dict(slopes=False)

        # Get section indices to run comparisons
        ind_st, ind_en, x, x = get_section(section, self.nROIs, self.roi_lr)

        # Calculate the meg connectivity data for each oscillation band
        for key in self.meg_con.keys():
            self.meg_con[key] = _mat_mult(self.meg_ROI_maps[key][ind_st:ind_en])

        # Initialize a dictionary to store data
        meg_stats = _init_meg_map_dict(length=2, slopes=False)

        # Get nROIs used in comparison
        nROIs_used = ind_en - ind_st

        # Extract anat range to use
        anat_comp = self.anat_con[ind_st:ind_en, ind_st:ind_en]

        # Calculate the correlations between each oscillation and anat data
        for key in meg_stats.keys():
            meg_stats[key][0], meg_stats[key][1] = sps.pearsonr(
                self.meg_con[key][np.triu_indices(nROIs_used, 1)],
                anat_comp[np.triu_indices(nROIs_used, 1)])

        # Attach the meg stats dictionary to object
        self.meg_stats = meg_stats

        # Print out results, if asked for
        if print_out:
            print('Anatomical data used is: ', self.anat_type)
            print('Correlation between MEG and anatomical connectivity: \n')

            # Loop through each oscillation, and print out R-val and p-val
            for key in self.meg_stats.keys():
                print(key)
                print('    R value: ', format(self.meg_stats[key][0], '1.4'))
                print('    P value: ', format(self.meg_stats[key][1], '1.4'))


#################################################################################################
############################ OMEGAMAPPIN - CL_MC - PRIVATE FUNCTIONS ############################
#################################################################################################


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
    names : list(str)
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


def _init_meg_map_dict(length=0, slopes=True):
    """Initialize a dictionary to store meg data.

    Parameters
    ----------
    length : int, optional
        If non-zero, length of zeros array to initialize.
            Defaults to 0, which initializes an empty array.
    slopes : boolean, optional
        Whether to include the 'Slopes' key in the dictonary.
            Defaults to True

    Returns
    -------
    meg_map : dictionary
        Dictionary with fields for MEG oscillation data.
    """

    # If length 0, initialize dict with empty arrays
    if length == 0:
        meg_map = dict([('Theta',     np.array([])),
                        ('Alpha',     np.array([])),
                        ('Beta',      np.array([])),
                        ('LowGamma',  np.array([])),
                        ('Slopes',    np.array([]))])

    # If given a number, initialize zeros arrays of given length
    else:
        meg_map = dict([('Theta',     np.zeros(length)),
                        ('Alpha',     np.zeros(length)),
                        ('Beta',      np.zeros(length)),
                        ('LowGamma',  np.zeros(length)),
                        ('Slopes',    np.zeros(length))])

    # If requested, drop the Slopes field
    if not slopes:
        meg_map.pop('Slopes')

    return meg_map


def _init_stat_dict():
    """Initialize a dictionary to store stat data for inter-data correlations.

    Returns
    -------
    out : dict
        A dictionary to store stats for terms/genes across all meg dat types.
    """

    out = dict([('TermsTheta',  np.array([])), ('TermsAlpha',    np.array([])),
                ('TermsBeta',   np.array([])), ('TermsLowGamma', np.array([])),
                ('GenesTheta',  np.array([])), ('GenesAlpha',    np.array([])),
                ('GenesBeta',   np.array([])), ('GenesLowGamma', np.array([])),
                ('TermsSlopes', np.array([])), ('GenesSlopes',   np.array([]))])

    return out


def _mat_mult(vec):
    """Multiply a vector element by element with itself.

    Parameters
    ----------
    vec : 1d array
        A vector to be multiplied by itself.

    Returns
    -------
    out : 2d array
        A matrix of the input vector multiplied by itself.
    """

    # Get length of first vector
    l = len(vec)

    # Initialize a matrix
    out = np.zeros([l, l])

    # Loop through vector, multiplying each element
    for i in range(0, l):
        for j in range(0, l):
            out[i, j] = vec[i] * vec[j]

    return out


def _make_list(dat_df):
    """Turns a dataframe into a list.

    Parameters
    ----------
    dat_df : DataFrame
        Input data, to be converted into a list.
    """

    # Get size of the data to
    [n_vert, n_dat] = dat_df.shape

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
    dat_in : list(tuple)
        A list of correlation results. Each tuple is (R-val, p-val).
    """

    # Check length of data
    n = len(dat_in)

    # Initializ vectors
    out_1 = np.zeros([n, 1])
    out_2 = np.zeros([n, 1])

    # Loop through and pull out data
    for i in range(n):
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

    # Return results
    return (corr_vals, p_vals)
