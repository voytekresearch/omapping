"""   """

import numpy as np
from py.test import raises

from om.core.io import *
from om.tests.utils import TestDB as TDB
from om.meg.single import MegData

##########################################################################################
##########################################################################################
##########################################################################################

def test_load_meg_psds():
    """   """

    tdb = TDB()
    pass

def test_save_foof_pickle():
    """   """

    tdb = TDB()

    foof_dat = [(1., np.array([5., 10.]), np.array([1., 2.]), np.array([1., 1.])),
                (1.5, np.array([10., 15.]), np.array([1., 2.]), np.array([1., 1.]))]

    save_foof_pickle(foof_dat, tdb.foof_path, 999)

    assert os.path.exists(os.path.join(tdb.foof_path, 'pickle', '999_Foof_Vertex.p'))

def test_load_foof_pickle():
    """   """

    tdb = TDB()

    results = load_foof_pickle(tdb.foof_path, 999)

    assert results

def test_save_foof_csv():
    """   """

    tdb = TDB()

    foof_dat = [(1., np.array([5., 10.]), np.array([1., 2.]), np.array([1., 1.])),
                (1.5, np.array([10., 15.]), np.array([1., 2.]), np.array([1., 1.]))]

    save_foof_csv(foof_dat, tdb.foof_path, 999)

    assert os.path.exists(os.path.join(tdb.foof_path, 'csv', '999_Slopes.csv'))
    assert os.path.exists(os.path.join(tdb.foof_path, 'csv', '999_Oscs.csv'))

def test_save_obj_pickle():

    tdb = TDB()

    dat = MegData(tdb, '')

    save_obj_pickle(dat, 'meg', 'test', db=tdb)

    assert True

    with raises(UnknownDataTypeError):
        save_obj_pickle(dat, 'BAD', 'test', db=tdb)

def test_load_obj_pickle():

    tdb = TDB()

    assert load_obj_pickle('meg', 'test', db=tdb)

    with raises(UnknownDataTypeError):
        load_obj_pickle('BAD', 'test', db=tdb)
