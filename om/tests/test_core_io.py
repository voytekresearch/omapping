"""Test for OM - core IO functions."""

import numpy as np
from py.test import raises

from om.core.io import *
from om.core.osc import Osc
from om.tests.utils import TestDB as TDB
from om.meg.single import MegSubj

##########################################################################################
##########################################################################################
##########################################################################################

def test_load_meg_psds():
    """   """

    tdb = TDB()
    pass

def test_save_fooof_pickle():
    """   """

    tdb = TDB()

    fooof_dat = [(1., np.array([5., 10.]), np.array([1., 2.]), np.array([1., 1.])),
                (1.5, np.array([10., 15.]), np.array([1., 2.]), np.array([1., 1.]))]

    save_fooof_pickle(fooof_dat, tdb.fooof_path, 999)

    assert os.path.exists(os.path.join(tdb.fooof_path, 'pickle', '999_fooof_Vertex.p'))

def test_load_fooof_pickle():
    """   """

    tdb = TDB()

    results = load_fooof_pickle(tdb.fooof_path, 999)

    assert results

def test_save_fooof_csv():
    """   """

    tdb = TDB()

    fooof_dat = [(1., np.array([5., 10.]), np.array([1., 2.]), np.array([1., 1.])),
                (1.5, np.array([10., 15.]), np.array([1., 2.]), np.array([1., 1.]))]

    save_fooof_csv(fooof_dat, tdb.fooof_path, 999)

    assert os.path.exists(os.path.join(tdb.fooof_path, 'csv', '999_Slopes.csv'))
    assert os.path.exists(os.path.join(tdb.fooof_path, 'csv', '999_Oscs.csv'))

def test_save_obj_pickle():

    tdb = TDB()

    dat = MegSubj(tdb, '')

    save_obj_pickle(dat, 'meg', 'test', db=tdb)

    assert True

    with raises(UnknownDataTypeError):
        save_obj_pickle(dat, 'BAD', 'test', db=tdb)

def test_load_obj_pickle():

    tdb = TDB()

    assert load_obj_pickle('meg', 'test', db=tdb)

    with raises(UnknownDataTypeError):
        load_obj_pickle('BAD', 'test', db=tdb)

def test_load_meg_list():

    tdb = TDB()
    osc = Osc(default=True)

    dat = load_meg_list(['test_v5', 'test_v5'], osc_bands_vert=True,
                        all_oscs=True, osc=osc, db=tdb, dat_source='')

    assert dat
