from __future__ import print_function, division

import numpy as np
from py.test import raises

from helper_test_funcs import load_test_meg_subj
from helper_test_funcs import TestDB as TDB

import om.cl.md_sing as md
from om.gen import OMDB, Osc

######################################################################################
##################### TESTS - OMEGAMAPPIN - CL_MD_SING - CLASSES #####################
######################################################################################

def test_meg_data():
    """   """

    db = OMDB()

    assert md.MegData(db, '')

###########################################################################################
################### TESTS - OMEGAMAPPIN - CL_MD_SING - PUBLIC FUNCTIONS ###################
###########################################################################################

def test_print_corrs_mat():

    rs_dat = np.array([[0.5, 0.5], [0.5, 0.5]])
    ps_dat = np.array([[0.1, 0.1], [0.1, 0.1]])
    labels = ['a', 'b']

    md.print_corrs_mat(rs_dat, ps_dat, labels)

    assert True

def test_print_corrs_vec():

    rs_dat = np.array([0.5, 0.5])
    ps_dat = np.array([0.1, 0.1])

    labels = ['a', 'b']
    desc = 'dat'

    md.print_corrs_vec(rs_dat, ps_dat, labels, desc)

    assert True

def test_save_md_pickle():

    tdb = TDB()

    dat = md.MegData(tdb, '')

    md.save_md_pickle(dat, 'test.p')

    assert True

def test_load_md_pickle():

    pass
    #dat = md.load_md_pickle('test.p')
    #assert True

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

    out_dat = md._get_single_osc(cens, pows, bws, osc_l, osc_h)

    assert np.array_equal(out_dat, np.array([9, 2, 1, 2]))

def test_get_single_osc_power_none():
    """   """

    cens = np.array([])
    pows = np.array([])
    bws = np.array([])

    out_cen, out_pow, out_bw = md._get_single_osc_power(cens, pows, bws)

    assert out_cen == 0.
    assert out_pow == 0.
    assert out_bw == 0.

def test_get_single_osc_power_one():

    cens = np.array([10.])
    pows = np.array([1.])
    bws = np.array([1.])

    out_cen, out_pow, out_bw = md._get_single_osc_power(cens, pows, bws)

    assert out_cen == 10.
    assert out_pow == 1.
    assert out_bw == 1.

def test_get_single_osc_power_mult():

    cens = np.array([10., 11.])
    pows = np.array([1., 2.])
    bws = np.array([1., 1.])

    out_cen, out_pow, out_bw = md._get_single_osc_power(cens, pows, bws)

    assert out_cen == 11.
    assert out_pow == 2.
    assert out_bw == 1.

def test_get_demo_csv_omega():
    """   """

    db = OMDB()

    subnum = 220216
    dat_source = 'OMEGA'

    sex, age = md._get_demo_csv(subnum, db.meg_path, dat_source)

    assert sex == 'M'
    assert age == 30

def test_get_demo_csv_hcp():
    """   """

    db = OMDB()

    subnum = 146129
    dat_source = 'HCP'

    sex, age = md._get_demo_csv(subnum, db.meg_path, dat_source)

    assert sex == 'M'
    assert age == 23.5

def test_osc_peak():
    """   """

    # Initialize data
    centers = np.array([1, 2, 3, 4])
    osc_low = 1.5
    osc_high = 3.5

    # Run function
    test_ans = md._osc_peak(centers, osc_low, osc_high)

    # Check the returned answer is close to the actual answer
    actual_ans = (2 + 3) / 2
    assert np.isclose(test_ans, actual_ans)

##########################################################################################
################### TESTS - OMEGAMAPPIN - CL_MD_SING - CLASS FUNCTIONS ###################
##########################################################################################

def test_cl_set_bands():

    tdb = TDB()

    meg_dat = md.MegData(tdb, '')

    osc = Osc(default=True)

    meg_dat.set_bands(osc)

def test_cl_set_bands_error():

    tdb = TDB()

    osc = Osc(default=True)

    meg_dat = md.MegData(tdb, '', osc=osc)

    meg_dat.bands_vertex = True

    with raises(md.InconsistentDataError):
        meg_dat.set_bands(osc)

def test_cl_import_foof():

    tdb = TDB()

    meg_dat = md.MegData(tdb, '')

    meg_dat.import_foof('test_v2', get_demo=False, load_type='pickle')

    exp_osc_count = np.array([2, 3])

    assert meg_dat.subnum == 'test_v2'
    assert meg_dat.n_psds == len(meg_dat.slopes) == 2
    assert meg_dat.centers.shape == meg_dat.powers.shape == meg_dat.bws.shape == (2, 8)
    assert np.array_equal(meg_dat.osc_count, exp_osc_count)
    assert meg_dat.has_data

def test_cl_import_foof_error():

    tdb = TDB()

    meg_dat = md.MegData(tdb, '')
    meg_dat.import_foof('test_v2', get_demo=False, load_type='pickle')

    with raises(md.InconsistentDataError):
        meg_dat.import_foof('test_v5', get_demo=False, load_type='pickle')

def test_cl_osc_bands_vertex():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.osc_bands_vertex()

    # TODO: ADD MORE PROPER CHECKS HERE!!!
    assert meg_dat.bands_vertex

def test_cl_all_oscs():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.all_oscs(verbose=True)

    # TODO: ADD MORE PROPER CHECKS HERE!!!
    assert meg_dat.all_osc

def test_cl_all_oscs_nan():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.bws[0, 0] = np.nan

    meg_dat.all_oscs(verbose=True)

    # TODO: ADD MORE PROPER CHECKS HERE!!!
    assert meg_dat.all_osc

def test_cl_peak_freq():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.all_oscs()

    meg_dat.peak_freq()

    # TODO: ADD MORE PROPER CHECKS HERE!!! CONSIDER HOW TO USE OSCS FOR PEAKS.
    assert meg_dat.peaks

def test_cl_peak_freq_errors():

    tdb = TDB()

    meg_dat = md.MegData(tdb, '')

    with raises(md.DataNotComputedError):
        meg_dat.peak_freq()

    meg_dat.all_osc = True

    with raises(md.DataNotComputedError):
        meg_dat.peak_freq()

def test_cl_calc_osc_param_corrs():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.all_oscs()

    a, b, c = meg_dat.calc_osc_param_corrs()

    # TODO: ADD MORE PROPER CHECKS HERE!!! CONSIDER HOW TO USE OSCS FOR PEAKS.
    assert a.any()

def test_cl_set_foof_viz():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.osc_bands_vertex()

    meg_dat.set_foof_viz()

    assert True