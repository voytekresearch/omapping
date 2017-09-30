"""Tests for OM - MEG decoding."""

from py.test import raises

from om.meg.decoding import *

##########################################################################################
##########################################################################################

def test_load_subjs():
    pass

def test_split_inds():
    """   """

    train_inds, test_inds = split_inds(100, 80, 20)

    assert len(set(train_inds)) == 80
    assert len(set(test_inds)) == 20

    # Check there is no overlap in indices
    assert set(train_inds).isdisjoint(set(test_inds))

    # Check error catching
    with raises(InconsistentDataError):
        assert split_inds(80, 80, 20)

def test_check_accuracy_single():
    pass

def test_check_accuracy_all():
    pass

def test_knn():
    pass