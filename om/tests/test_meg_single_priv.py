"""   """

import numpy as np

import om.meg.single as single
from om.core.db import OMDB

##########################################################################################
################## TESTS - OMEGAMAPPIN - CL_MD_SING - PRIVATE FUNCTIONS ##################
##########################################################################################

def test_get_single_osc():
    """   """

    cens = np.array([3, 6, 9, 12, 15])
    pows = np.array([1, 1, 2, 1, 1])
    bws = np.array([1, 1, 1, 1, 1])

    osc_l = 5
    osc_h = 10

    out_dat = single._get_single_osc(cens, pows, bws, osc_l, osc_h)

    assert np.array_equal(out_dat, np.array([9, 2, 1, 2]))

def test_get_single_osc_power_none():
    """   """

    cens = np.array([])
    pows = np.array([])
    bws = np.array([])

    out_cen, out_pow, out_bw = single._get_single_osc_power(cens, pows, bws)

    assert out_cen == 0.
    assert out_pow == 0.
    assert out_bw == 0.

def test_get_single_osc_power_one():

    cens = np.array([10.])
    pows = np.array([1.])
    bws = np.array([1.])

    out_cen, out_pow, out_bw = single._get_single_osc_power(cens, pows, bws)

    assert out_cen == 10.
    assert out_pow == 1.
    assert out_bw == 1.

def test_get_single_osc_power_mult():

    cens = np.array([10., 11.])
    pows = np.array([1., 2.])
    bws = np.array([1., 1.])

    out_cen, out_pow, out_bw = single._get_single_osc_power(cens, pows, bws)

    assert out_cen == 11.
    assert out_pow == 2.
    assert out_bw == 1.

def test_get_demo_csv_omega():
    """   """

    db = OMDB()

    subnum = 220216
    dat_source = 'OMEGA'

    sex, age = single._get_demo_csv(subnum, db.meg_path, dat_source)

    assert sex == 'M'
    assert age == 30

def test_get_demo_csv_hcp_unres():
    """   """

    db = OMDB()

    subnum = 146129
    dat_source = 'HCP'

    sex, age = single._get_demo_csv(subnum, db.meg_path, dat_source, use_restricted=False)

    assert sex == 'M'
    assert age == 23.5

def test_osc_peak_all():
    """   """

    centers = np.array([6, 9, 10, 12, 15])
    osc_low = 8
    osc_high = 13

    test_ans = single._osc_peak_all(centers, osc_low, osc_high)

    assert np.isclose(test_ans, 10.33, rtol=0.01)

    test_ans = single._osc_peak_all(centers, osc_low, osc_high, 'median')
    assert np.isclose(test_ans, 10.00, rtol=0.01)
