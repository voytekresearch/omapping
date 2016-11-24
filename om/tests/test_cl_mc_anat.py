from __future__ import print_function, division

import numpy as np
from py.test import raises

from om.gen import OMDB
import om.cl.mc_anat as mc
from helper_test_funcs import TestDB as TDB

##############################################################################
################# TESTS - OMEGAMAPPIN - CL_MC_ANAT - CLASSES #################
##############################################################################

def test_roi():
    """   """

    roi = mc.ROI()

    assert roi

def test_roi_set_labels():
    """   """

    roi = mc.ROI()

    labels = ['a', 'b', 'c']
    lrs = ['l', 'r', 'l']

    roi.set_labels(labels, lrs)

    assert roi.labels
    assert roi.lrs
    assert roi.n_rois

def test_roi_set_labels_error():
    """   """

    roi = mc.ROI()

    labels = ['a', 'b', 'c']
    lrs = ['l', 'r']

    with raises(mc.InconsistentDataError):
        assert roi.set_labels(labels, lrs)

def test_roi_set_verts():
    """   """

    roi = mc.ROI()

    verts = np.array([1, 2, 3])
    roi.set_verts(verts)

    assert np.any(roi.verts)

def test_roi_set_comment():
    """   """

    roi = mc.ROI()

    comment = 'test'
    roi.set_comment(comment)

    assert roi.comment

def test_check_consistency1():
    """Test check consistency, with label and lr data."""

    roi = mc.ROI()

    labels = ['a', 'b', 'c']
    lrs = ['l', 'r', 'l']

    roi.set_labels(labels, lrs)

    roi.check_consistency()

    roi.labels.pop()

    with raises(mc.InconsistentDataError):
        assert roi.check_consistency()

def test_check_consistency2():
    """Test check consistency, with label, lr and verts data."""

    roi = mc.ROI()

    labels = ['a', 'b', 'c']
    lrs = ['l', 'r', 'l']

    verts = np.array([[1, 2], [3, 4], [5, 6]])

    roi.set_labels(labels, lrs)
    roi.set_verts(verts)

    roi.check_consistency()

    verts = np.array([[1, 2], [3, 4]])
    roi.set_verts(verts)

    with raises(mc.InconsistentDataError):
        roi.check_consistency()

def test_mc_anat():
    """   """

    db = OMDB()

    assert mc.MapCompAnat(db)

####################################################################################
################### TESTS - OMEGAMAPPIN - CL_MC_ANAT - FUNCTIONS ###################
####################################################################################

def test_extract_lr_anat():
    """   """

    anat_labels = ['left_caudal_anterior_cingulate',
                   'right_caudal_middle_frontal']

    labels, lrs = mc._extract_lr(anat_labels, 'anat')

    assert len(labels) == 2
    assert len(lrs) == 2
    assert labels[0] == 'caudalanteriorcingulate'
    assert labels[1] == 'caudalmiddlefrontal'
    assert lrs[0] is 'L'
    assert lrs[1] is 'R'

def test_extract_lr_elec():
    """   """

    elec_labels = ['parsopercularis R', 'postcentral L',]

    labels, lrs = mc._extract_lr(elec_labels, 'elec')

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

    label_1 = mc._clean_label(label_1, lr_1, 'anat')
    label_2 = mc._clean_label(label_2, lr_2, 'anat')

    assert label_1 == 'caudalanteriorcingulate'
    assert label_2 == 'caudalanteriorcingulate'

def test_clean_label_elec():
    """   """

    label = 'parsopercularis R'
    lr = 'R'

    label = mc._clean_label(label, lr, 'elec')

    assert label == 'parsopercularis'

def test_mat_mult():
    """   """

    dat = np.array([1, 2, 3])

    res = mc._mat_mult(dat)
    exp = np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]])

    assert res.shape == (3, 3)
    assert np.all(res == exp)

####################################################################################
################# TESTS - OMEGAMAPPIN - CL_MC_ANAT - ROI FUNCTIONS #################
####################################################################################

def test_load_anat_maps():

    tdb = TDB()

    map_comp = mc.MapCompAnat(tdb)

    map_comp.load_anat_maps('test_anat.mat', 'test')

    assert map_comp.anat.loaded

def test_load_elec_rois():

    tdb = TDB()

    map_comp = mc.MapCompAnat(tdb)

    map_comp.load_elec_rois('test_scout.mat')

    assert map_comp.elec.loaded

def test_load_elec_rois_default():

    db = OMDB()

    map_comp = mc.MapCompAnat(db)
    map_comp.load_elec_rois()

    assert map_comp.elec.loaded

def test_align_rois():

    tdb = TDB()

    map_comp = mc.MapCompAnat(tdb)

    map_comp.load_anat_maps('test_anat.mat', 'test')
    map_comp.load_elec_rois('test_scout.mat')

    map_comp.align_rois()

    assert map_comp.rois.loaded

def test_align_rois_errors():

    tdb = TDB()

    map_comp = mc.MapCompAnat(tdb)

    with raises(mc.DataNotComputedError):
        assert map_comp.align_rois()

def test_conv_meg_rois():

    tdb = TDB()

    map_comp = mc.MapCompAnat(tdb)

    map_comp.load_meg_maps('test_meg')

    map_comp.load_anat_maps('test_anat.mat', 'test')
    map_comp.load_elec_rois('test_scout.mat')

    map_comp.align_rois()

    map_comp.conv_meg_rois()

    assert True

def test_comp_meg_anat():

    tdb = TDB()

    map_comp = mc.MapCompAnat(tdb)

    map_comp.load_meg_maps('test_meg')

    map_comp.load_anat_maps('test_anat.mat', 'test')
    map_comp.load_elec_rois('test_scout.mat')

    map_comp.align_rois()

    map_comp.conv_meg_rois()

    map_comp.comp_meg_anat()

    assert True
