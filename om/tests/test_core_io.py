
import numpy as np

from om.tests.utils import TestDB as TDB

from om.core.io import *

#import om.cl.md_sing as md
from om.meg.single import MegData

##
##
##

def test_load_meg_psds():
    """   """

    tdb = TDB()
    pass

def test_load_foof_pickle():
    pass

def test_save_foof_pickle():
    """   """

    tdb = TDB()

    foof_dat = [(1., np.array([5., 10.]), np.array([1., 2.]), np.array([1., 1.])),
                (1.5, np.array([10., 15.]), np.array([1., 2.]), np.array([1., 1.]))]

    save_foof_pickle(foof_dat, tdb.foof_path, 999)

    assert os.path.exists(os.path.join(tdb.foof_path, 'pickle', '999_Foof_Vertex.p'))

def test_save_foof_csv():

    tdb = TDB()

    foof_dat = [(1., np.array([5., 10.]), np.array([1., 2.]), np.array([1., 1.])),
                (1.5, np.array([10., 15.]), np.array([1., 2.]), np.array([1., 1.]))]

    save_foof_csv(foof_dat, tdb.foof_path, 999)

    assert os.path.exists(os.path.join(tdb.foof_path, 'csv', '999_Slopes.csv'))
    assert os.path.exists(os.path.join(tdb.foof_path, 'csv', '999_Oscs.csv'))

def test_save_md_pickle():

    tdb = TDB()

    dat = MegData(tdb, '')

    save_md_pickle(dat, 'test', db=tdb)

    assert True

def test_load_md_pickle():

    tdb = TDB()

    assert load_md_pickle('test', db=tdb)
