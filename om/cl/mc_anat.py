"""MODULE DOCSTRING - TO FILL IN

"""

from __future__ import print_function, division
import os
import numpy as np
import scipy.io as sio
import scipy.stats.stats as sps

# Import custom om code
from om.gen import get_section, DataNotComputedError
from om.cl.mc_base import MapCompBase, _init_meg_map_dict

################################################################
########## OMEGAMAPPIN - MAP COMPARE - ANAT - CLASSES ##########
################################################################

class MapCompAnat(MapCompBase):
    """Class for storing and comparing spatial topographies in rois."""

    def __init__(self, db):
        """

        Parameters
        ----------
        db : ?
            xx
        """

        # Inherit from MapComp() class
        MapCompBase.__init__(self, db)

        # Initialize var to store number of rois
        self.n_rois = int()

        # Add vars to save roi data from anat data
        self.anat_roi_names = list()
        self.anat_roi_lrs = list()
        self.anat_n_rois = int()
        self.anat_con = np.ndarray(0)
        self.anat_type = str()

        # Add vars to save roi data from MEG data
        self.elec_roi_names = list()
        self.elec_roi_lrs = list()
        self.elec_n_rois = int
        self.elec_roi_verts = list()

        # Initialize list to store roi labels
        self.roi_labels = list()
        self.roi_verts = list()
        self.roi_lr = list()

        # Initialize var to store MEG roi data
        self.meg_roi_maps = dict()

        # Initialize var to store MEG connectivity data
        self.meg_con = dict()

        # Initialize var to store meg stats data
        self.meg_stats = dict()

        # Add boolean for whether anat data is loaded
        self.anat_loaded = False
        self.elec_loaded = False
        self.rois_aligned = False
        self.meg_rois = False


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

        # Get full path for the anat mat file
        anat_mat_file = os.path.join(self.db.maps_anat_path, anat_file_name)

        # Load the anat data
        dat = sio.loadmat(anat_mat_file)

        # Pull out data from mat dictionary
        roi_names = dat['roi_labels'].tolist()
        self.anat_con = dat['connectivity']

        # Get number of rois
        self.anat_n_rois = len(roi_names)

        # Loop through and fix roi labels
        for roi in range(0, self.anat_n_rois):
            self.anat_roi_names.append(str(roi_names[roi][0][0]))

        # Update which type of anatomy data is loaded
        self.anat_type = anat_type

        # Update boolean that anat data is loaded
        self.anat_loaded = True


    def load_elec_rois(self, roi_file_name=None):
        """Load the roi file for the MEG data.

        Parameters
        ----------
        self : MapComproi() object
            Object for storing and comparing map data, in roi format.
        roi_file_name : str, optional
            File name (including path) to the matlab scout definition file.
        """

        # Set default scouts file, if an alternative was not provided
        if not roi_file_name:
            roi_file_name = 'scout_Desikan-Killiany_68.mat'

        # Load roi dat from mat file
        dat = sio.loadmat(os.path.join(self.db.maps_scouts_path, roi_file_name),
                          struct_as_record=True)

        # Pull out data from mat file
        scouts = np.squeeze(dat['Scouts'])

        # Check how many rois there are
        self.elec_n_rois = len(scouts)

        # Initiliaze vars to store scout names and vertices
        sc_names = list()
        sc_verts = list()

        # Loop through, pull out names and verts into lists
        for i in range(0, self.elec_n_rois):
            scout = scouts[i]
            sc_verts.append(scout[0])
            sc_names.append(str(scout[3]))

            # Drop brackets in scout name
            sc_names[i] = sc_names[i][3:-2]

        # Attach data to object
        self.elec_roi_names = sc_names
        self.elec_roi_verts = sc_verts

        # Update boolean that elec data is loaded
        self.elec_loaded = True


    def align_rois(self):
        """Align rois used and names between anat and meg rois.

        NOTES
        -----
        xxxxxxxx
        """

        # Check if rois loaded - return if not
        if (not self.elec_roi_names) or (not self.anat_roi_names):
            raise DataNotComputedError('One or Both rois not loaded! Cant proceed!')

        # Elec L/Rs - standardize names
        for roi in range(0, self.elec_n_rois):

            # Check if left side scout
            if self.elec_roi_names[roi][-1] is 'L':
                self.elec_roi_lrs.append('L')

            # Check if right side scout
            elif self.elec_roi_names[roi][-1] is 'R':
                self.elec_roi_lrs.append('R')

            else:
                pass

            # Drop the L/R from the names
            self.elec_roi_names[roi] = self.elec_roi_names[roi][:-2]

        # Anat L/Rs - standardize names
        for roi in range(0, self.anat_n_rois):

            # Check if left side scout
            if self.anat_roi_names[roi][0] is 'l':
                self.anat_roi_lrs.append('L')
                self.anat_roi_names[roi] = self.anat_roi_names[roi][5:]

            # Check if right side scout
            elif self.anat_roi_names[roi][0] is 'r':
                self.anat_roi_lrs.append('R')
                self.anat_roi_names[roi] = self.anat_roi_names[roi][6:]

            else:
                pass

            # Drop the underscores
            self.anat_roi_names[roi] = self.anat_roi_names[roi].replace("_", "")

        # Loop through and line up scout names
        for i in range(0, self.anat_n_rois):

            # Grab current roi from anat roi list
            cur_roi = self.anat_roi_names[i]
            cur_lr = self.anat_roi_lrs[i]

            # Loop through elec rois to match up
            for j in range(0, self.elec_n_rois):

                # Check if current elec roi matches current roi
                if self.elec_roi_names[j] == cur_roi:

                    # Check if L/R is right
                    if self.elec_roi_lrs[j] == cur_lr:

                        # Same side - add to overall list
                        self.roi_labels.append(cur_roi)
                        self.roi_lr.append(cur_lr)

                        # Add vertices to overall roi list
                        self.roi_verts.append(self.elec_roi_verts[j])

        # Check how many rois there are
        self.n_rois = len(self.roi_labels)

        # Set boolean that rois have been aligned
        self.rois_aligned = True


    def conv_meg_rois(self):
        """Convert MEG data to rois.

        NOTES
        -----
        This XXXXX
        """

        # Initialize dict for current roi data
        roi_meg_dat = _init_meg_map_dict(self.bands.keys(), self.n_rois)

        # Loop through all rois
        for roi in range(0, self.n_rois):

            # Add current roi data to dict
            # Loop through all oscs
            for key in self.meg_maps.keys():

                #
                cur_verts = np.squeeze(self.roi_verts[roi] - 1)
                n_verts = len(cur_verts)

                #
                temp_dat = self.meg_maps[key][cur_verts]

                #
                roi_meg_dat[key][roi] = (sum(temp_dat) / n_verts)

        # Add the current roi data to object
        self.meg_roi_maps = roi_meg_dat

        # Update boolean that meg data has been converted to rois
        self.meg_rois = True


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
        ind_st, ind_en, _, _ = get_section(section, self.n_rois, self.roi_lr)

        # Calculate the meg connectivity data for each oscillation band
        for key in self.meg_con.keys():
            self.meg_con[key] = _mat_mult(self.meg_roi_maps[key][ind_st:ind_en])

        # Initialize a dictionary to store data
        meg_stats = _init_meg_map_dict(self.bands.keys(), length=2)

        # Get n_rois used in comparison
        n_rois_used = ind_en - ind_st

        # Extract anat range to use
        anat_comp = self.anat_con[ind_st:ind_en, ind_st:ind_en]

        # Calculate the correlations between each oscillation and anat data
        for key in meg_stats.keys():
            meg_stats[key][0], meg_stats[key][1] = sps.pearsonr(
                self.meg_con[key][np.triu_indices(n_rois_used, 1)],
                anat_comp[np.triu_indices(n_rois_used, 1)])

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


################################################################################################
######################## OMEGAMAPPIN - CL_MC_ANAT - FUNCTIONS (PRIVATE) ########################
################################################################################################

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
