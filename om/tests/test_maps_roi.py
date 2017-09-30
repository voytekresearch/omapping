"""Tests for OM - maps roi."""

import numpy as np
from py.test import raises

from om.maps.roi import *

##############################################################################
###################### TESTS - OMEGAMAPPIN - MAPS - ROI ######################
##############################################################################

def test_roi():
    """   """

    roi = ROI()

    assert roi

def test_roi_set_labels():
    """   """

    roi = ROI()

    labels = ['a', 'b', 'c']
    lrs = ['l', 'r', 'l']

    roi.set_labels(labels, lrs)

    assert roi.labels
    assert roi.lrs
    assert roi.n_rois

def test_roi_set_labels_error():
    """   """

    roi = ROI()

    labels = ['a', 'b', 'c']
    lrs = ['l', 'r']

    with raises(InconsistentDataError):
        assert roi.set_labels(labels, lrs)

def test_roi_set_verts():
    """   """

    roi = ROI()

    verts = np.array([1, 2, 3])
    roi.set_verts(verts)

    assert np.any(roi.verts)

def test_roi_set_comment():
    """   """

    roi = ROI()

    comment = 'test'
    roi.set_comment(comment)

    assert roi.comment

def test_check_consistency1():
    """Test check consistency, with label and lr data."""

    roi = ROI()

    labels = ['a', 'b', 'c']
    lrs = ['l', 'r', 'l']

    roi.set_labels(labels, lrs)

    roi.check_consistency()

    roi.labels.pop()

    with raises(InconsistentDataError):
        assert roi.check_consistency()

def test_check_consistency2():
    """Test check consistency, with label, lr and verts data."""

    roi = ROI()

    labels = ['a', 'b', 'c']
    lrs = ['l', 'r', 'l']

    verts = np.array([[1, 2], [3, 4], [5, 6]])

    roi.set_labels(labels, lrs)
    roi.set_verts(verts)

    roi.check_consistency()

    verts = np.array([[1, 2], [3, 4]])
    roi.set_verts(verts)

    with raises(InconsistentDataError):
        roi.check_consistency()
