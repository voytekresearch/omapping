"""Tests for OM - meg twin."""

from om.core.osc import Osc
from om.tests.utils import TestDB as TDB
from om.tests.utils import load_test_meg_pair

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

def test_comp_peak_freq():
    """Tests the comp_peak_freq() function.

    Note: this calls (and tests) sub-function _comp_pf_pair().
    """

    dat = load_test_meg_pair(osc_bands_vert=True)

    out = comp_peak_freq([dat, dat])

    assert out

def test_comp_osc_space():
    """Tests the comp_osc_space() function.

    Note: this calls (and tests) sub-function _comp_space_pair().
    """

    dat = load_test_meg_pair(osc_bands_vert=True)

    out = comp_osc_space([dat, dat])

    assert out

def test_comp_osc_param():
    """Tests the comp_osc_param() function.

    Note: this calls (and tests) sub-function _comp_osc_pair().
    """

    dat = load_test_meg_pair(osc_bands_vert=True)

    out = comp_osc_param([dat, dat], 0)

    assert out

def test_comp_slope():
    """Tests the comp_slope() function.

    Note: this calls (and tests) sub-function _comp_sl_pair().
    """

    dat = load_test_meg_pair()

    out = comp_slope([dat, dat])

    assert out

def test_print_funcs():
    """Tests the functions to print out results run."""

    vec_dat = np.array([1, 2])
    corr_dat = np.array([[1, 1], [2, 2]])
    labels = ['aa', 'bb']

    print_twin_results_vec(vec_dat, labels)
    print_twin_results_corr(corr_dat, labels)

    assert True
