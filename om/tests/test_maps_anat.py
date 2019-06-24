"""Tests for OM - maps anat"""

from py.test import raises

import numpy as np

from om.maps.anat import *
from om.core.db import OMDB
from om.tests.utils import TestDB as TDB
from om.tests.utils import load_test_anat

from om.maps.anat import _extract_lr, _clean_label, _mat_mult

###################################################################################################
###################################################################################################

def test_mc_anat():
    """   """

    db = OMDB()

    assert MapCompAnat(db)

def test_load_anat_maps():

    tdb = TDB()

    map_comp = MapCompAnat(tdb)

    map_comp.load_anat_maps('test_anat.mat', 'test')

    assert map_comp.anat.loaded

def test_load_elec_rois():

    tdb = TDB()

    map_comp = MapCompAnat(tdb)

    map_comp.load_elec_rois('test_scout.mat')

    assert map_comp.elec.loaded

def test_load_elec_rois_default():

    db = OMDB()

    map_comp = MapCompAnat(db)
    map_comp.load_elec_rois()

    assert map_comp.elec.loaded

def test_align_rois():

    map_comp = load_test_anat(load_anat=True, load_scout=True)

    map_comp.align_rois()

    assert map_comp.rois.loaded

def test_align_rois_errors():

    tdb = TDB()

    map_comp = MapCompAnat(tdb)

    with raises(DataNotComputedError):
        assert map_comp.align_rois()

def test_conv_meg_rois():


    map_comp = load_test_anat(load_meg=True, load_scout=True, load_anat=True,
                              align=True)

    map_comp.conv_meg_rois()

    assert True

def test_calc_meg_con():
    """   """

    map_comp = load_test_anat(load_meg=True, load_scout=True, load_anat=True,
                              align=True, convert=True)

    map_comp.calc_meg_con()

    assert True

def test_comp_meg_anat():

    map_comp = load_test_anat(load_meg=True, load_scout=True, load_anat=True,
                              align=True, convert=True, calc_meg=True)

    map_comp.comp_meg_anat(print_out=False)

    assert True

def test_check_comp():

    map_comp = load_test_anat()

    with raises(DataNotComputedError):
        map_comp.check_comp()

    map_comp = load_test_anat(load_meg=True, load_scout=True, load_anat=True,
                          align=True, convert=True, calc_meg=True)

    map_comp.comp_meg_anat(print_out=False)

    map_comp.check_comp()

    map_comp.comp_meg_anat(print_out=True)

    assert True

######################################################################################
################### TESTS - OMEGAMAPPIN - ANAT - PRIVATE FUNCTIONS ###################
######################################################################################


def test_extract_lr_anat():
    """   """

    anat_labels = ['left_caudal_anterior_cingulate',
                   'right_caudal_middle_frontal']

    labels, lrs = _extract_lr(anat_labels, 'anat')

    assert len(labels) == 2
    assert len(lrs) == 2
    assert labels[0] == 'caudalanteriorcingulate'
    assert labels[1] == 'caudalmiddlefrontal'
    assert lrs[0] is 'L'
    assert lrs[1] is 'R'

def test_extract_lr_elec():
    """   """

    elec_labels = ['parsopercularis R', 'postcentral L',]

    labels, lrs = _extract_lr(elec_labels, 'elec')

    assert len(labels) == 2
    assert len(lrs) == 2
    assert labels[0] == 'parsopercularis'
    assert labels[1] == 'postcentral'
    assert lrs[0] is 'R'
    assert lrs[1] is 'L'

def test_clean_label_anat():
    """   """

    label_1 = 'left_caudal_anterior_cingulate'
    label_2 = 'right_caudal_anterior_cingulate'
    lr_1 = 'L'
    lr_2 = 'R'

    label_1 = _clean_label(label_1, lr_1, 'anat')
    label_2 = _clean_label(label_2, lr_2, 'anat')

    assert label_1 == 'caudalanteriorcingulate'
    assert label_2 == 'caudalanteriorcingulate'

def test_clean_label_elec():
    """   """

    label = 'parsopercularis R'
    lr = 'R'

    label = _clean_label(label, lr, 'elec')

    assert label == 'parsopercularis'

def test_mat_mult():
    """   """

    dat = np.array([1, 2, 3])

    res = _mat_mult(dat)
    exp = np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]])

    assert res.shape == (3, 3)
    assert np.all(res == exp)
