"""Oscillation band defintion object for OM.

Notes
-----
 - The Osc object uses the OrderedDict object, which is a special dictionary object
which enforces and ensures a consistent order of items within the dictionary.
More information here: https://docs.python.org/2/library/collections.html#collections.OrderedDict
"""
from __future__ import print_function

from types import StringType
from collections import OrderedDict
import numpy as np

from om.core.errors import InconsistentDataError

######################################################################################
############################## OMEGAMAPPIN - CORE - OSC ##############################
######################################################################################

class Osc(object):
    """Class to hold definitions of oscillation bands.

    Attributes
    ----------
    bands : dict
        Dictionary of oscillation band definitions.
    labels : list of str
        Labels for all oscillation bands.
    n_bands : int
        The number of oscillation bands that are defined.
    """

    def __init__(self, default=False, input_bands=None):
        """Initialize the Osc() object.

        Parameters
        ----------
        default : boolean, optional (default = False)
            Whether to use the default oscillation bands.
        input_bands : dict, optional (default = None)
            A dictionary of oscillation bands to use.

        Notes:
        - If supplied, an input_bands will over-ride the default bands option,
            even if it is set as True.
        ^This behaviour is unhelpful. Fix up.
        """

        # Initialize bands as a dictionary
        self.bands = OrderedDict()

        # Add running total of how many bands are defined
        self.n_bands = int()

        # If requested use the default oscillation bands
        if default:
            self.add_band('Theta', (3, 8))
            self.add_band('Alpha', (8, 13))
            self.add_band('Beta', (13, 30))
            self.add_band('LowGamma', (30, 40))

        # If supplied, use the given dictionary of oscillation bands
        if input_bands:

            # Print a warning if adding custom bands on top of the default ones
            if default:
                print('WARNING: Input bands will be added to the default bands.')

            # Initialize lists
            keys = []
            lower_bounds = []

            # Loop through the provided dictionary
            for key in input_bands:

                # Check that provided bands are legal
                _check_band(key, input_bands[key])

                # Keep a running list of the band names and lower bound
                keys.append(key)
                lower_bounds.append(input_bands[key][0])

            # Sort the bands, and add them in ascending order
            sort_inds = np.argsort(lower_bounds)
            for ind in sort_inds:
                self.add_band(keys[ind], input_bands[keys[ind]])


    def add_band(self, band_name, band_lims):
        """Add a new oscillation band definition.

        Parameters
        ----------
        band_name : str
            The name of the new oscillation band.
        band_lims : tuple(float, float)
            The lower and upper frequency limit of the band.
        """

        # Check that band definition is properly formatted
        _check_band(band_name, band_lims)

        # Add the given band to oscillation definition
        self.bands[band_name] = band_lims
        self.n_bands += 1

        # Update labels
        self.labels = self.bands.keys()


    def rm_band(self, rm_band):
        """Remove a previously defined oscillation band.

        Parameters
        ----------
        rm_band : str
            Band name to remove from oscillation band definitions.
        """

        # Remove requested band from oscillation definition
        self.bands.pop(rm_band)
        self.n_bands -= 1

##################################################################################################
##################################################################################################
##################################################################################################

def check_bands(osc_lst):
    """Check that a list of oscillation band definitions are all the same. If so, return bands.

    Parameters
    ----------
    osc_lst : dict
        Oscillation band definitions to compare. Should be dict from Osc.bands.

    Returns
    -------
    OrderedDict
        Oscillation band oscillations.
    """

    # Check that all oscillation definitions provided are the same
    if not all(x == osc_lst[0] for x in osc_lst):
        raise InconsistentDataError('Oscillation definitions are inconsistent.')

    return osc_lst[0]


def CheckBands(func):
    """Decorator function to check that all oscillation band definitions are consistent.

    Notes
    -----
    - This decorator requires that 'dat', a variable with a list of MegData objects be
        the first argument for any function it is wrapped around.
    - It also requires the wrapped function to take an argument 'bands', but it need not
        be given, as this decorator will provide it.
    """

    def wrapper(dat, *args, **kwargs):
        """Wrapper that checks band definitions, and passes through given inputs & bands."""

        bands = check_bands([subj.bands for subj in dat])

        return func(dat, bands=bands, *args, **kwargs)

    return wrapper

###################################################################################################
###################################################################################################
###################################################################################################

def _check_band(band_name, band_limits):
    """Check that a proposed band definition is properly formatted.

    Parameters
    ----------
    band_name : str
        The name of the new oscillation band.
    band_limits : tuple(float, float)
        The lower and upper frequency limit of the band.

    Raises
    ------
    InconsistentDataError
        If oscillation band definition is not properly formatted.
    """

    # Check that band name is a string
    if not isinstance(band_name, StringType):
        raise InconsistentDataError('Band name definition is not a string.')

    # Check that band limits has the right size
    if not len(band_limits) == 2:
        raise InconsistentDataError('Band limit definition is not the right size.')

    # Safety check that limits are in correct order
    if not band_limits[0] < band_limits[1]:
        raise InconsistentDataError('Band limits are incorrect.')
