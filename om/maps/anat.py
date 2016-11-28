"""MODULE DOCSTRING - TO FILL IN

"""

from __future__ import print_function, division

import os
import numpy as np
import scipy.io as sio
import scipy.stats.stats as sps

# Import custom om code
from om.maps.base import MapCompBase, _init_meg_map_dict
from om.maps.roi import ROI
from om.core.utils import get_section
from om.core.errors import DataNotComputedError

####################################################################################
#################### OMEGAMAPPIN - MAP COMPARE - ANAT - CLASSES ####################
####################################################################################

class MapCompAnat(MapCompBase):
    """Class for storing and comparing spatial topographies in ROIs.

    This class inherits from MapCompBase() and inherits all it's attributes.

    Attributes
    ----------
    rois : ROI() object
        ROI data aligned across anat & elec.
    anat : ROI() object
        ROI information from anatomical data.
    elec : ROI() object
        ROI information from MEG data
    anat_con : 2d array
        Matrix of anatomical connectivity data.
    meg_con : dict
        MEG physiological matrices for all oscillation bands.
    meg_roi_maps : dict
        MEG data collapsed into ROIs. Average data, in a given band, across the ROI.
    stats : dict
        Correlation results from comparing anatomy to electrophysiology.
    """

    def __init__(self, db):
        """

        Parameters
        ----------
        db : OMDB() object
            Database object for omegamappin project.
        """

        # Inherit from MapComp() class
        MapCompBase.__init__(self, db)


        # Initialize vars to store ROI data
        self.rois = ROI()
        self.anat = ROI()
        self.elec = ROI()

        # ?
        self.anat_con = np.ndarray(0)

        # Initialize var to store MEG roi data
        self.meg_roi_maps = dict()

        # Initialize var to store MEG connectivity data
        self.meg_con = dict()

        # Initialize var to store correlation results
        self.stats = dict()


    def load_anat_maps(self, anat_file_name, anat_type):
        """Load the spatial maps of anatomilcal data.

        Parameters
        ----------
        self : MapComproi() object
            Object for storing and comparing map data, in roi format.
        anat_file_name : str
            File name of anat data file.
        anat_type : str
            Indicates which type of anat data is loaded.
        """

        # Initialize ROI object to store anat ROI data
        anat = ROI()

        # Set up file for anatomy data, and load
        anat_mat_file = os.path.join(self.db.maps_path, 'Anat', anat_file_name)
        dat = sio.loadmat(anat_mat_file)

        # Pull out data from mat dictionary
        label_dat = dat['roi_labels'].tolist()
        con_dat = dat['connectivity']

        # Pull out label data
        extracted_labels = []
        for label in label_dat:
            extracted_labels.append(str(label[0][0]))

        # Pull out L/R info from labels
        labels, lrs = _extract_lr(extracted_labels, 'anat')

        # Add data to anat ROI object
        anat.set_labels(labels, lrs)
        anat.set_comment('anat_' + anat_type)

        # Attach anat ROI object, and anat data to MapComp object
        self.anat = anat
        self.anat_con = con_dat


    def load_elec_rois(self, roi_file_name=None):
        """Load the roi file for the MEG data.

        Parameters
        ----------
        roi_file_name : str, optional
            File name (including path) to the matlab scout definition file.
        """

        # Initialize ROI object to store anat ROI data
        elec = ROI()

        # Set default scouts file, if an alternative was not provided
        if not roi_file_name:
            roi_file_name = 'scout_Desikan-Killiany_68.mat'

        # Load roi dat from mat file
        dat = sio.loadmat(os.path.join(self.db.maps_path, 'Scouts', roi_file_name),
                          struct_as_record=True)

        # Pull out data from mat file
        scouts = np.squeeze(dat['Scouts'])

        # Initiliaze vars to store scout names and vertices
        sc_names = list()
        sc_verts = list()

        # Loop through, pull out names and verts into lists
        for ind, scout in enumerate(scouts):

            sc_verts.append(scout[0])
            sc_names.append(str(scout[3]))

            # Drop brackets in scout name
            sc_names[ind] = sc_names[ind][3:-2]

        # Extract L/R data from names
        labels, lrs = _extract_lr(sc_names, 'elec')

        # Set label and vertex data
        elec.set_labels(labels, lrs)
        elec.set_verts(sc_verts)
        elec.set_comment('elec')

        # Attach anat ROI object to MapComp object
        self.elec = elec


    def align_rois(self):
        """Align rois used and names between anat and meg rois.

        NOTES
        -----
        xxxxxxxx
        """

        # Check if rois loaded - return if not
        if (not self.elec.loaded) or (not self.anat.loaded):
            raise DataNotComputedError('One or Both rois not loaded! Cant proceed!')

        # Initialize roi object for aligned object
        rois = ROI()

        # Loop through and line up scout names
        # Loops through anat ROIs, finding corresponding ROI in elec
        for anat_ind, anat_label in enumerate(self.anat.labels):

            anat_lr = self.anat.lrs[anat_ind]

            for elec_ind, elec_label in enumerate(self.elec.labels):

                elec_lr = self.elec.lrs[elec_ind]

                # Check labels match
                if anat_label == elec_label:

                    # Check if L/R is right
                    if anat_lr == elec_lr:

                        #
                        rois.labels.append(anat_label)
                        rois.lrs.append(anat_lr)

                        #
                        rois.verts.append(self.elec.verts[elec_ind])

        #
        rois.n_rois = len(rois.labels)
        rois.loaded = True

        #
        rois.check_consistency()

        # Add
        self.rois = rois


    def conv_meg_rois(self):
        """Convert MEG data to rois.

        NOTES
        -----
        This XXXXX
        """

        # Initialize dict for current roi data
        roi_meg_dat = _init_meg_map_dict(self.bands.keys(), self.rois.n_rois)

        # Loop through all rois
        for ind, cur_verts in enumerate(self.rois.verts):

            # Add current roi data to dict
            # Loop through all oscs
            for band in self.meg_maps.keys():

                n_verts = len(cur_verts)

                temp_dat = self.meg_maps[band][cur_verts]

                roi_meg_dat[band][ind] = (sum(temp_dat) / n_verts)

        # Add the current roi data to object
        self.meg_roi_maps = roi_meg_dat


    def comp_meg_anat(self, section='all', print_out=True):
        """Compare anatomical connectivity to oscillation data.

        Parameters
        ----------
        self : MapComproi() object
            Object for storing and comparing map data, in roi format.
        section : {'all, 'left', 'right'}
            Which section of data to compare.
        print_out : boolean, optional (default = True)
            Whether to print out the stats results.
        """

        # Initialize the dictionary to store MEG connectivity data
        self.meg_con = _init_meg_map_dict(self.bands.keys())

        # Get section indices to run comparisons
        ind_st, ind_en, _, _ = get_section(section, self.rois.n_rois, self.rois.lrs)

        # Calculate the meg connectivity data for each oscillation band
        for key in self.meg_con.keys():
            self.meg_con[key] = _mat_mult(self.meg_roi_maps[key][ind_st:ind_en])

        # Initialize a dictionary to store data
        stats = _init_meg_map_dict(self.bands.keys(), length=2)

        # Get n_rois used in comparison
        n_rois_used = ind_en - ind_st

        # Extract anat range to use
        anat_comp = self.anat_con[ind_st:ind_en, ind_st:ind_en]

        # Calculate the correlations between each oscillation and anat data
        for key in stats.keys():
            stats[key][0], stats[key][1] = sps.pearsonr(
                self.meg_con[key][np.triu_indices(n_rois_used, 1)],
                anat_comp[np.triu_indices(n_rois_used, 1)])

        # Attach the stats dictionary to object
        self.stats = stats

        # Print out results, if asked for
        if print_out:
            print('Anatomical data used is: ', self.anat.comment)
            print('Correlation between MEG and anatomical connectivity: \n')

            # Loop through each oscillation, and print out R-val and p-val
            for key in self.stats.keys():
                print(key)
                print('    R value: ', format(self.stats[key][0], '1.4'))
                print('    P value: ', format(self.stats[key][1], '1.4'))


############################################################################################
###################### OMEGAMAPPIN - CL_MC_ANAT - FUNCTIONS (PRIVATE) ######################
############################################################################################

def _extract_lr(labels_in, dat_type):
    """

    Parameters
    ----------
    labels : list of str
        xx
    dat_type : {'anat', 'elec'}
        xx

    Returns
    -------
    labels : list of str
        xx
    lrs : list of str
        xx
    """

    # Initialize lists to collect data
    labels_out = []
    lrs = []

    # Set index to check, based on data type
    if dat_type is 'anat':
        ch_ind = 0
    elif dat_type is 'elec':
        ch_ind = -1

    #
    for label in labels_in:

        #
        if label[ch_ind].lower() == 'l':
            lr = 'L'
        elif label[ch_ind].lower() == 'r':
            lr = 'R'

        #
        label = _clean_label(label, lr, dat_type)

        #
        labels_out.append(label)
        lrs.append(lr)

    return labels_out, lrs


def _clean_label(label, lr, dat_type):
    """

    Parameters
    ----------
    label : str
        xx
    lr : str
        xx
    dat_type : {'anat', 'elec'}
        xx
    """

    #
    if dat_type is 'elec':
        label = label[:-2]

    #
    if dat_type is 'anat':
        if lr is 'L':
            label = label[5:]
        elif lr is 'R':
            label = label[6:]

        # For anat data, drop underscores in the name
        label = label.replace("_", "")

    return label


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
    vec_len = len(vec)

    # Initialize a matrix
    out = np.zeros([vec_len, vec_len])

    # Loop through vector, multiplying each element
    for i in range(0, vec_len):
        for j in range(0, vec_len):
            out[i, j] = vec[i] * vec[j]

    return out
