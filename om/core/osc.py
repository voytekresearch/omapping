"""   """

from om.core.errors import InconsistentDataError

####################################################################################
############################# OMEGAMAPPIN - CORE - OSC #############################
####################################################################################

class Osc(object):
    """Class to hold definition of oscillation bands.

    Attributes
    ----------
    bands : dict
        Dictionary of oscillation band definitions.
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
        """

        # Initialize bands as a dictionary
        self.bands = dict()

        # If requested use the default oscillation bands
        if default:
            self.bands = dict({'Theta': (3, 8),
                               'Alpha': (8, 13),
                               'Beta': (13, 30),
                               'LowGamma': (30, 40)})

        # If supplied, use the given dictionary of oscillation bands
        if input_bands:
            self.bands = input_bands


    def add_band(self, band_name, band_lims):
        """Add a new oscillation band definition.

        Parameters
        ----------
        band_name : str
            The name of the new oscillation band.
        band_lims : tuple(float, float)
            The lower and upper frequency limit of the band.

        Raises
        ------
        InconsistentDataError
            If oscillation band limits given do not work.
        """

        # Safety check that limits are in correct order
        if not band_lims[0] < band_lims[1]:
            raise InconsistentDataError('Band limits are incorrect.')

        # Add the given band to oscillation definition
        self.bands[band_name] = band_lims


    def rm_band(self, old_band):
        """Remove a previously defined oscillation band.

        Parameters
        ----------
        old_band : str
            Band name to remove from oscillation band definitions.
        """

        # Remove requested band from oscillation definition
        self.bands.pop(old_band)
