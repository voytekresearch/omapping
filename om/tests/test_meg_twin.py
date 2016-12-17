"""DOCSTRING"""

from om.core.osc import Osc
from om.tests.utils import TestDB as TDB

from om.meg.twin import *

##################################################################################
##################################################################################
##################################################################################

def test_get_twin_data():
    """   """

    a, b, c, d = get_twin_data()

    assert True

def test_match_twins():
    """   """

    twin_dat = np.array([[121, 111, 333], [122, 222, 444], [123, 111, 333]])

    twins_out, single_out = match_twins(twin_dat)

    assert twins_out[0] == (121, 123)
    assert single_out[0] == tuple([122])

def test_check_complete_pairs():
    """   """

    twin_pairs = [(111, 114), (112, 116), (110, 117)]
    av_dat = [111, 112, 113, 114, 115]

    complete_pairs = check_complete_pairs(twin_pairs, av_dat)

    assert len(complete_pairs) == 1
    assert complete_pairs[0] == (111, 114)

def test_rm_twin_pairs():
    """   """

    all_pairs = [(111, 112), (111, 113), (112, 113), (113, 114)]
    twin_pairs = [(111, 113), (113, 114)]

    non_twins = rm_twin_pairs(all_pairs, twin_pairs)

    assert len(non_twins) == 2
    assert set([(111, 112) , (112, 113)]) == set(non_twins)

def test_compare_pair():
    """   """

    tdb = TDB()
    osc = Osc(default=True)

    corr_dat = compare_pair(['test_v5', 'test_v5'], osc=osc, db=tdb, dat_source='')

    assert np.any(corr_dat)

def test_compare_slope():
    """   """

    tdb = TDB()

    corr_dat = compare_slope(['test_v5', 'test_v5'], db=tdb, dat_source='')

    assert np.any(corr_dat)

def test_print_twin_results():
    """   """

    corr_dat = np.array([[1, 1], [2, 2]])
    labels = ['aa', 'bb']

    print_twin_results(corr_dat, labels)

    assert True
