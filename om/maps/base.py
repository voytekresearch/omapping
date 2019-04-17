"""MODULE DOCSTRING - TO FILL IN"""

import os
import pickle
import numpy as np

#######################################################################################
############################## OMEGAMAPPIN - MAPS - BASE ##############################
#######################################################################################

class MapCompBase(object):
    """Class for storing and comparing spatial topographies."""

    def __init__(self, db):
        """Initialize object with omegamappin database.

        Parameters
        ----------
        db : OMDB() object
            Database object for omegamappin project.
        """

        # Add database object
        self.db = db

        # Initialize a dictionary to store maps of meg data (oscillation bands)
        self.meg_maps = dict()
        self.bands = dict()

        # Initialize a dictionary to store exponent map
        self.exponent_map = dict({'Exponents': np.array([])})

        # Initialize booleans that keep track of what is loaded
        self.oscs_loaded = False
        self.exponents_loaded = False


    def load_meg_maps(self, osc_file):
        """Load the spatial maps of MEG data (oscillation bands).

        Parameters
        ----------
        osc_file : str
            File name of the pickle file with oscillations data.
        """

        # Get the full path for the file name
        osc_maps_file = os.path.join(self.db.maps_path, 'Oscs', osc_file + '.p')

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


    def load_exponent_map(self, exponent_file):
        """Load the spatial map of MEG exponent data.

        Parameters
        ----------
        exponent_file : str
            File name of the pickle file with exponent data.
        """

        # Get the full path for the file name
        exponents_map_file = os.path.join(self.db.maps_path, 'Exponents', exponent_file + '.p')

        # Load data from pickle file
        dat_in = pickle.load(open(exponents_map_file, 'rb'))

        # Pull out the exponent data
        self.exponent_map['Exponents'] = dat_in['exponents']

        # Update boolean that exponents are loaded
        self.exponents_loaded = True


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
    meg_map : dict
        Dictionary with fields for MEG oscillation data.
    """

    # Initialize dictionary
    meg_map = dict()

    # Add oscillation bands
    for band in bands:
        meg_map[band] = np.zeros(length)

    return meg_map
