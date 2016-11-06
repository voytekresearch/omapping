from __future__ import print_function, division
import os
import pickle
import numpy as np
import scipy.io as sio

# Import custom om code
from om.gen import *

############################################################################################
############################ OMEGAMAPPIN - CL_MC_BASE - CLASSES ############################
############################################################################################

class MapCompBase(object):
    """Class for storing and comparing spatial topographies."""

    def __init__(self, db):
        """

        Parameters
        ----------
        db : xx
            xx
        """

        # Add database object
        self.db = db

        # Initialize a dictionary to store maps of meg data (oscillation bands)
        self.meg_maps = dict()
        self.bands = dict()

        # Initialize a dictionary to store slope map
        self.slope_map = dict({'Slopes': np.array([])})

        # Initialize booleans that keep track of what is loaded
        self.oscs_loaded = False
        self.slopes_loaded = False


    def load_meg_maps(self, osc_file):
        """Load the spatial maps of MEG data (oscillation bands).

        Parameters
        ----------
        self : MapComp() object
            Object for storing and comparing map data.
        osc_file : str, optional
            File name of the pickle file with oscillations data.
        """

        # Get the full path for the file name
        osc_maps_file = os.path.join(self.db.maps_oscs_path, osc_file + '.p')

        # Load data from pickle file
        dat_in = pickle.load(open(osc_maps_file, 'rb'))

        # Get the oscillation bands used in current maps
        self.bands = dat_in['bands']

        # Initialize the var to store meg map data
        self.meg_maps = _init_meg_map_dict(self.bands.keys())

        # Pull out oscillation band data
        for band in self.bands:
            self.meg_maps[band] = dat_in['osc_dat'][band]

        # Update boolean that oscs are loaded
        self.oscs_loaded = True


    def load_slope_map(self, slope_file):
        """Load the spatial map of MEG slope data.

        Parameters
        ----------
        slope_file : str
            File name of the pickle file with slope data.
        """

        # Get the full path for the file name
        slopes_map_file = os.path.join(self.db.maps_slopes_path, slope_file + '.p')

        # Load data from pickle file
        dat_in = pickle.load(open(slopes_map_file, 'rb'))

        # Pull out the slope data
        self.slope_map['Slopes'] = dat_in['slopes']

        # Update boolean that slopes are loaded
        self.slopes_loaded = True


################################################################################################
######################## OMEGAMAPPIN - CL_MC_BASE - FUNCTIONS (PRIVATE) ########################
################################################################################################

def _init_meg_map_dict(bands, length=0):
    """Initialize a dictionary to store meg data.

    Parameters
    ----------
    bands : list of str
        Oscillation bands to initialize.
    length : int, optional (default = 0)
        If non-zero, length of zeros array to initialize.

    Returns
    -------
    meg_map : dictionary
        Dictionary with fields for MEG oscillation data.
    """

    # Initialize dictionary
    meg_map = dict()

    # Add oscillation bands
    for band in bands:
        meg_map[band] = np.zeros(length)

    return meg_map
