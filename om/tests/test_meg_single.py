"""   """

import numpy as np
from py.test import raises

from om.meg.single import *
from om.meg.single import _get_single_osc, _get_single_osc_power, _get_demo_csv, _osc_peak_all

from om.core.db import OMDB
from om.core.osc import Osc
from om.tests.utils import TestDB as TDB
from om.tests.utils import load_test_meg_subj

# TODO: Update tests to use updated boolean structure of available data

########################################################################################
#################### TESTS - OMEGAMAPPIN - MEG - SINGLE - FUNCTIONS ####################
########################################################################################

def test_print_corrs_mat():

    rs_dat = np.array([[0.5, 0.5], [0.5, 0.5]])
    ps_dat = np.array([[0.1, 0.1], [0.1, 0.1]])
    labels = ['a', 'b']

    print_corrs_mat(rs_dat, ps_dat, labels)

    assert True

def test_print_corrs_vec():

    rs_dat = np.array([0.5, 0.5])
    ps_dat = np.array([0.1, 0.1])

    labels = ['a', 'b']
    desc = 'dat'

    print_corrs_vec(rs_dat, ps_dat, labels, desc)

    assert True

######################################################################################
######################### TESTS - OMEGAMAPPIN - MEG - SINGLE #########################
######################################################################################

def test_meg_data():
    """   """

    db = OMDB()

    assert MegData(db, '')

def test_set_bands():

    tdb = TDB()

    meg_dat = MegData(tdb, '')

    osc = Osc(default=True)

    meg_dat.set_bands(osc)

def test_set_bands_error():

    tdb = TDB()

    osc = Osc(default=True)

    meg_dat = MegData(tdb, '', osc=osc)

    meg_dat.has_vertex_bands = True

    with raises(InconsistentDataError):
        meg_dat.set_bands(osc)

def test_import_foof():

    tdb = TDB()

    meg_dat = MegData(tdb, '')

    meg_dat.import_foof('test_v2', get_demo=False, load_type='pickle')

    exp_osc_count = np.array([2, 3])

    assert meg_dat.subnum == 'test_v2'
    assert meg_dat.n_psds == len(meg_dat.slopes) == 2
    assert meg_dat.centers.shape == meg_dat.powers.shape == meg_dat.bws.shape == (2, 8)
    assert np.array_equal(meg_dat.osc_count, exp_osc_count)
    assert meg_dat.has_data

def test_import_foof_error():

    tdb = TDB()

    meg_dat = MegData(tdb, '')
    meg_dat.import_foof('test_v2', get_demo=False, load_type='pickle')

    with raises(InconsistentDataError):
        meg_dat.import_foof('test_v5', get_demo=False, load_type='pickle')

def test_all_oscs():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.all_oscs(verbose=True)

    # TODO: ADD MORE PROPER CHECKS HERE!!!
    assert meg_dat.has_all_osc

def test_all_oscs_nan():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.bws[0, 0] = np.nan

    meg_dat.all_oscs(verbose=True)

    # TODO: ADD MORE PROPER CHECKS HERE!!!
    assert meg_dat.has_all_osc

def test_osc_bands_vertex():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.osc_bands_vertex()

    # TODO: ADD MORE PROPER CHECKS HERE!!!
    assert meg_dat.has_vertex_bands

def tests_osc_bands_vertex_error():

    tdb = TDB()

    meg_dat = MegData(tdb, '')

    with raises(DataNotComputedError):
        meg_dat.osc_bands_vertex()

def test_peak_freq_all():

    meg_dat = load_test_meg_subj('test_v2', all_oscs=True)

    meg_dat.peak_freq(dat='all')

    # TODO: ADD MORE PROPER CHECKS HERE!!! CONSIDER HOW TO USE OSCS FOR PEAKS.
    assert meg_dat.peaks

    # Check error
    meg_dat = load_test_meg_subj('test_v2')

    with raises(DataNotComputedError):
        meg_dat.peak_freq(dat='all')

def test_peak_freq_band():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.osc_bands_vertex()

    meg_dat.peak_freq(dat='band')

    assert meg_dat.peaks

    # Check median
    meg_dat.peak_freq(dat='band', avg='median')

    # Check error
    meg_dat = load_test_meg_subj('test_v2')

    with raises(DataNotComputedError):
        meg_dat.peak_freq(dat='band')

def test_peak_freq_errors():

    tdb = TDB()

    meg_dat = MegData(tdb, '')

    with raises(DataNotComputedError):
        meg_dat.peak_freq(dat='all')

    meg_dat.has_all_osc = True

    with raises(DataNotComputedError):
        meg_dat.peak_freq(dat='all')

def test_calc_osc_param_corrs():

    meg_dat = load_test_meg_subj('test_v2', all_oscs=True)

    a, b, c = meg_dat.calc_osc_param_corrs()

    # TODO: ADD MORE PROPER CHECKS HERE!!! CONSIDER HOW TO USE OSCS FOR PEAKS.
    assert a.any()

def test_calc_osc_param_corrs_error():

    meg_dat = load_test_meg_subj('test_v2')

    with raises(DataNotComputedError):
        meg_dat.calc_osc_param_corrs()

def test_set_foof_viz():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.osc_bands_vertex()

    meg_dat.set_foof_viz()

    assert True

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

    out_dat = _get_single_osc(cens, pows, bws, osc_l, osc_h)

    assert np.array_equal(out_dat, np.array([9, 2, 1, 2]))

def test_get_single_osc_power_none():
    """   """

    cens = np.array([])
    pows = np.array([])
    bws = np.array([])

    out_cen, out_pow, out_bw = _get_single_osc_power(cens, pows, bws)

    assert out_cen == 0.
    assert out_pow == 0.
    assert out_bw == 0.

def test_get_single_osc_power_one():

    cens = np.array([10.])
    pows = np.array([1.])
    bws = np.array([1.])

    out_cen, out_pow, out_bw = _get_single_osc_power(cens, pows, bws)

    assert out_cen == 10.
    assert out_pow == 1.
    assert out_bw == 1.

def test_get_single_osc_power_mult():

    cens = np.array([10., 11.])
    pows = np.array([1., 2.])
    bws = np.array([1., 1.])

    out_cen, out_pow, out_bw = _get_single_osc_power(cens, pows, bws)

    assert out_cen == 11.
    assert out_pow == 2.
    assert out_bw == 1.

def test_get_demo_csv_omega():
    """   """

    db = OMDB()

    subnum = 220216
    dat_source = 'OMEGA'

    sex, age = _get_demo_csv(subnum, db.meg_path, dat_source)

    assert sex == 'M'
    assert age == 30

def test_get_demo_csv_hcp_unres():
    """   """

    db = OMDB()

    subnum = 146129
    dat_source = 'HCP'

    sex, age = _get_demo_csv(subnum, db.meg_path, dat_source, use_restricted=False)

    assert sex == 'M'
    assert age == 23.5

def test_osc_peak_all():
    """   """

    centers = np.array([6, 9, 10, 12, 15])
    osc_low = 8
    osc_high = 13

    test_ans = _osc_peak_all(centers, osc_low, osc_high)

    assert np.isclose(test_ans, 10.33, rtol=0.01)

    test_ans = _osc_peak_all(centers, osc_low, osc_high, 'median')
    assert np.isclose(test_ans, 10.00, rtol=0.01)
