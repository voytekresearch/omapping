from __future__ import print_function, division
import os
import csv
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr
from om.gen import *

###################################################################################
######################## OMEGAMAPPIN - MAP COMPARE CLASSES ########################
###################################################################################

class MapComp():
    """Class for storing and comparing spatial topographies."""

    def __init__(self):

        # Set root path for where all map data is stored
        self.maps_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/Maps/'

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

        # Initialize a dictionary to store all the R-value results from spatial correlations
        self.corrs = dict([('TermsTheta',  np.array([])), ('TermsAlpha',    np.array([])),
                           ('TermsBeta',   np.array([])), ('TermsLowGamma', np.array([])),
                           ('GenesTheta',  np.array([])), ('GenesAlpha',    np.array([])),
                           ('GenesBeta',   np.array([])), ('GenesLowGamma', np.array([])),
                           ('TermsSlopes', np.array([])), ('GenesSlopes',   np.array([]))
                          ])

        # Initialize a dictionary to store all the p-value results from spatial correlations
        self.p_vals = dict([('TermsTheta',  np.array([])),  ('TermsAlpha',    np.array([])),
                            ('TermsBeta',   np.array([])),  ('TermsLowGamma', np.array([])),
                            ('GenesTheta',  np.array([])),  ('GenesAlpha',    np.array([])),
                            ('GenesBeta',   np.array([])),  ('GenesLowGamma', np.array([])),
                            ('TermsSlopes', np.array([])),  ('GenesSlopes',   np.array([]))
                            ])

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
        osc_files = clean_file_list(os.listdir(self.oscs_path), '.npz')
        slope_files = clean_file_list(os.listdir(self.slopes_path), '.npz')
        gene_files = clean_file_list(os.listdir(self.genes_path), '.csv')
        term_files = clean_file_list(os.listdir(self.terms_path), '.csv')

        # If asked for, print out lists of files
        if print_files:
            print('Oscillation Files: \n', '\n'.join(osc_files), '\n')
            print('Slope Files: \n', '\n'.join(slope_files), '\n')
            print('Terms Files: \n', '\n'.join(term_files), '\n')
            print('Genes Files: \n', '\n'.join(gene_files), '\n')

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


    def load_gene_maps(self, genes_file_names):
        """Load the spatial maps of gene data.

        Note:
            Input must be a list. If gene data in single file, use single item list.
            If list contains multiple files, these files will be loaded and
                concatenated to form the full gene maps.
            Order: First file should be first set of files.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data. 
        genes_file_names : list (str)
            list of files containing gene data
        """
        
        # If one file given, load this as the gene map
        if len(genes_file_names) == 1:
            self.gene_maps = pd.read_csv(genes_file_names[0], header=None)

        # If multiple files, load them all and concatenate
        else:
            # Loop through all files given
            for i in range(0, len(genes_file_names)):
                # If first file, initialize as first part of the gene map
                if i == 0:
                    genes_csv = os.path.join(self.genes_path, genes_file_names[0])
                    self.gene_maps = pd.read_csv(genes_csv, header=None)
                # For all subsequent files, concatenate to end of gene map
                else:
                    genes_csv = os.path.join(self.genes_path, genes_file_names[i])
                    temp_df = pd.read_csv(genes_csv, header=None)
                    self.gene_maps = pd.concat([self.gene_maps, temp_df])

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


    def calc_corrs(self, dat_type, meg_dat):
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
        """

        # Check with data type and set data accordingly
        if dat_type is 'Terms':
            n_comps = self.n_terms
            dat_df = self.term_maps
        elif dat_type is 'Genes':
            n_comps = self.n_genes
            dat_df = self.gene_maps
        else:
            print("Improper Data Type. Fix it.")
            return

        # Get the specified meg map
        try:
            meg_map = self.meg_maps[meg_dat]
        except KeyError:
            print('MEG Data not understood. Fix it.')

        # Initialize vectors to store correlation results
        corr_vals = np.zeros([n_comps, 1])
        p_vals = np.zeros([n_comps, 1])

        # Loop through all comparisons to run
        for comp in range(0, n_comps):

            # Pull out specific data (single term or gene)
            dat = np.array(dat_df.ix[:, comp])

            # Get inds of data that contains numbers
            inds_non_nan = [i for i in range(len(dat)) if np.isnan(dat[i]) == False]

            # Calculate correlation between data and meg map
            [corr_vals[comp], p_vals[comp]] = pearsonr(dat[inds_non_nan], meg_map[inds_non_nan])

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
        top :
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
            print("Data type not understood. Fix it.")

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

            # Update boolean that genes are not loaded
            self.genes_loaded = True

        # Otherwise, data type was not understood
        else:
            print('Data type not understood. Try again.')


    def plot_corrs(self, dat_type, meg_dat):
        """Plot the R and p values of specified correlation results.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data. 
        dat_type : str
            Data type to plot corrs for (Terms or Genes). 
        meg_dat : str
            Specific meg map to plot corrs for. 

        """

        # Check that asked for correlations have been computed
        #if not self.corrs[dat_type + meg_dat]:
        #    print("Those correlations not calculated. Quitting.")
        #    return

        # Get the specified data
        meg_corrs = np.squeeze(self.corrs[dat_type + meg_dat])
        meg_p_vals = np.squeeze(self.p_vals[dat_type + meg_dat])
        
        # Initialize subplot figure
        fig, ax = plt.subplots(1, 2, figsize=(12, 6))

        # Set supertitle for the plot
        plt.suptitle('Corr Stats for ' + dat_type + ' ' + meg_dat, fontsize=20, fontweight='bold')

        # Plot R values
        ax[0].plot(meg_corrs, '.')
        ax[0].set_title('R Values', {'fontsize': 16, 'fontweight': 'bold'})
        ax[0].set_ylim([-0.5, 0.5])

        # Plot the p values
        ax[1].plot(meg_p_vals, '.')
        ax[1].set_title('P Values', {'fontsize': 16, 'fontweight': 'bold'})


    def save_corrs(self, dat_type, meg_dat, save_as_npz=True, save_as_csv=True):
        """Save out the correlation results.

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data. 
        dat_type : str
            Data type to save corrs for ('Terms' or 'Genes'). 
        meg_dat : str
            MEG data to save corrs for. 
        save_as_npz : boolean, optional
            Whether to save an npz file. Default is True. 
        save_as_csv : boolean, optional
            Whether to save a csv file. Default is True. 
        """

        # Check which type of data and set names accordingly
        if dat_type is 'Terms':
            names = self.term_names
        elif dat_type is 'Genes':
            names = self.gene_names
        else:
            print("Data type not understood. Fix it.")

        # Check that asked for correlations have been computed
        #if not self.corrs[dat_type + meg_dat]:
        #    print("Those correlations not calculated. Quitting.")
        #    return

        # Get the correlation data of interest
        meg_corrs = np.squeeze(self.corrs[dat_type + meg_dat])
        meg_p_vals = np.squeeze(self.p_vals[dat_type + meg_dat])

        # Set basic file name to save data as
        file_name = 'Corrs_' + dat_type + meg_dat
        save_path = os.path.join('/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/Corrs/', dat_type)

        # Save a numpy npz file
        if save_as_npz:

            # Set up specific name/path for npz file
            npz_path = os.path.join(save_path, 'npz', '')
            outfile_npz = npz_path + file_name + '.npz'

            # Save out npz file
            np.savez(outfile_npz, dat_type=dat_type, meg_dat=meg_dat, names=names, 
                     meg_corrs=meg_corrs, meg_p_vals=meg_p_vals)

        # Save a csv file
        if save_as_csv:

            # Set up specific name/path for csv file
            csv_path = os.path.join(save_path, 'csv', '')
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

    def __init__(self):

        # Inherit from MapComp() class
        MapComp.__init__(self)

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
        roi_file_name : ?
            xx
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

        XXXX
        """

        # Check if ROIs loaded - return if not
        if (not self.elec_roi_names) or (not self.anat_roi_names):
            print('One or Both ROIs not loaded! Cant proceed!')
            return

        # Elec L/Rs - standardize names
        for r in range(0, self.elec_nROIs):

            # Check if left side scout
            if(self.elec_roi_names[r][-1] is 'L'):
                self.elec_roi_lrs.append('L')
            # Check if right side scout
            elif(self.elec_roi_names[r][-1] is 'R'):
                self.elec_roi_lrs.append('R')
            else:
                pass

            # Drop the L/R from the names
            self.elec_roi_names[r] = self.elec_roi_names[r][:-2]

        # Anat L/Rs - standardize names
        for r in range(0, self.anat_nROIs):

            # Check if left side scout
            if(self.anat_roi_names[r][0] is 'l'):
                self.anat_roi_lrs.append('L')
                self.anat_roi_names[r] = self.anat_roi_names[r][5:]
            # Check if right side scout
            elif(self.anat_roi_names[r][0] is 'r'):
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
                    if (self.elec_roi_lrs[j] == cur_lr):

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


    def comp_meg_anat(self, print_out=True):
        """Compare anatomical connectivity to oscillation data. 

        Parameters
        ----------
        print_out : boolean, optional
            Whether to print out the stats results.
                Defaults to True
        """

        # Initialize the dictionary to store MEG connectivity data
        self.meg_con = _init_meg_map_dict(slopes=False)

        # Calculate the meg connectivity data for each oscillation band
        for key in self.meg_con.keys():
            self.meg_con[key] = _mat_mult(self.meg_ROI_maps[key])

        # Initialize a dictionary to store data 
        meg_stats = _init_meg_map_dict(length=2, slopes=False)

        # Calculate the correlations between each oscillation and anat data
        for key in meg_stats.keys():
            meg_stats[key][0], meg_stats[key][1] = pearsonr(
                self.meg_con[key][np.triu_indices(self.nROIs, 1)],
                self.anat_con[np.triu_indices(self.nROIs, 1)])

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


    def plot_mat(self, dat, osc=None, section='all'):
        """Plot the connectivity matrix. 

        Parameters
        ----------
        self : MapCompROI() object
            Object for storing and comparing map data, in ROI format.
        dat : str
            Which data type to plot. 
                Options: 'anat', 'meg'
        osc : str, optional
            If data is MEG, which oscillation to plot. 
                Options: 'Theta', 'Alpha', 'Beta', 'LowGamma'. Default is None.
        section : str
            Which part of the matrix to plot. 
                Options: 'all', 'left', 'right', 'lr', 'rl'
        """
        
        # Initialize figure
        f = plt.figure(figsize=[12, 12])
        ax = f.add_subplot(111)

        # Set section to plot
        if section is 'all':
            ind_st_a = ind_st_b = 0
            ind_en_a = ind_en_b = self.nROIs
        elif section is 'left':
            ind_st_a = ind_st_b = self.roi_lr.index('L')
            ind_en_a = ind_en_b = _find_last(self.roi_lr, 'L')
        elif section is 'right':
            ind_st_a = ind_st_b = self.roi_lr.index('R')
            ind_en_a = ind_en_b =  _find_last(self.roi_lr, 'R')
        elif section is 'lr':
            ind_st_a = self.roi_lr.index('L')
            ind_en_a = _find_last(self.roi_lr, 'L')
            ind_st_b = self.roi_lr.index('R')
            ind_en_b = _find_last(self.roi_lr, 'R')
        elif section is 'rl':
            ind_st_a = self.roi_lr.index('R')
            ind_en_a = _find_last(self.roi_lr, 'R')
            ind_st_b = self.roi_lr.index('L')
            ind_en_b = _find_last(self.roi_lr, 'L') 
        else:
            print('Section range unclear!')
            return

        # Create the figure
        if dat is 'anat':
            #m = ax.matshow(self.anat_con)
            m = ax.imshow(self.anat_con[ind_st_a:ind_en_a, ind_st_b:ind_en_b], interpolation='none')
        elif dat is 'meg':
            #m = ax.matshow(self.meg_con[osc])
            m = ax.imshow(self.meg_con[osc][ind_st_a:ind_en_a, ind_st_b:ind_en_b], interpolation='none')

        # Add colour bar as a legend
        f.colorbar(m)

        # Add a title to the plot
        if dat is 'anat':
            plt.title('Anatomy Connectivity', fontsize=22, fontweight='bold')
        elif dat is 'meg':
            plt.title('MEG Connectivity', fontsize=22, fontweight='bold')

        # Display the figure
        plt.show()


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
    """Initialize a dictionary to store meg data 

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
