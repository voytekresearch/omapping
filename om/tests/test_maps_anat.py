"""   """

import numpy as np
from py.test import raises

#from om.gen import OMDB
from om.core.db import OMDB

#import om.cl.mc_anat as mc
from om.maps.anat import *

#from helper_test_funcs import TestDB as TDB
from om.tests.utils import TestDB as TDB

#####################################################################################
######################### TESTS - OMEGAMAPPIN - MAPS - ANAT #########################
#####################################################################################

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

    tdb = TDB()

    map_comp = MapCompAnat(tdb)

    map_comp.load_anat_maps('test_anat.mat', 'test')
    map_comp.load_elec_rois('test_scout.mat')

    map_comp.align_rois()

    assert map_comp.rois.loaded

def test_align_rois_errors():

    tdb = TDB()

    map_comp = MapCompAnat(tdb)

    with raises(DataNotComputedError):
        assert map_comp.align_rois()

def test_conv_meg_rois():

    tdb = TDB()

    map_comp = MapCompAnat(tdb)

    map_comp.load_meg_maps('test_meg')

    map_comp.load_anat_maps('test_anat.mat', 'test')
    map_comp.load_elec_rois('test_scout.mat')

    map_comp.align_rois()

    map_comp.conv_meg_rois()

    assert True

def test_comp_meg_anat():

    tdb = TDB()

    map_comp = MapCompAnat(tdb)

    map_comp.load_meg_maps('test_meg')

    map_comp.load_anat_maps('test_anat.mat', 'test')
    map_comp.load_elec_rois('test_scout.mat')

    map_comp.align_rois()

    map_comp.conv_meg_rois()

    map_comp.comp_meg_anat()

    assert True
