"""DOCSTRING"""

import numpy as np

from om.core.errors import InconsistentDataError

###################################################################################################
###################################################################################################

class ROI(object):
    """Class to hold data about ROIs.

    Attributes
    ----------
    n_rois : int
        The number of ROIs.
    labels : list of str
        Label for each ROI.
    verts : list of array
        List of vertices included in the ROI.
    lr : list of str
        L/R for each ROI.
    comment : str
        Comment about the ROIs.
    """

    def __init__(self):
        """Initialize ROI object."""

        # Initialize fields to store ROI data
        self.n_rois = int()
        self.labels = list()
        self.lrs = list()
        self.verts = list()
        self.comment = str()
        self.loaded = bool()


    def set_labels(self, labels, lrs):
        """Set scout labels for current ROI object.

        Parameters
        ----------
        labels : ?
            xx
        lrs : ?
            xx
        """

        # Check labels and lrs are the same size
        if len(labels) != len(lrs):
            raise InconsistentDataError('Input labels and lrs do not match!')

        # Attach data to object
        self.labels = labels
        self.lrs = lrs
        self.n_rois = len(labels)
        self.loaded = True


    def set_verts(self, verts):
        """Set the vertices that define the scouts for current ROI object.

        Parameters
        ----------
        verts : ?
            xx
        """

        for vert in verts:
            self.verts.append(np.squeeze(vert - 1))


    def set_comment(self, comment):
        """Add a comment to current ROI object.

        Parameters
        ----------
        comment : str
            A comment regarding the current data.
        """

        self.comment = comment


    def check_consistency(self):
        """Check that ROI data is internally consistent."""

        # Check that labels and L/R definitions have the same number of elements
        if not self.n_rois == len(self.labels) == len(self.lrs):
            raise InconsistentDataError('Discrepancy in number of labels and/or L/Rs.')

        # If the ROI has vertex data defined, check it is consistent with labels
        if self.verts:
            if not self.n_rois == len(self.verts):
                raise InconsistentDataError('Discrepancy in the number of vertices.')
