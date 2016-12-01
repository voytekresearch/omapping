"""   """

import numpy as np
from py.test import raises

from om.meg.single import *
from om.core.db import OMDB
from om.core.osc import Osc
from om.tests.utils import TestDB as TDB
from om.tests.utils import load_test_meg_subj

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

def test_cl_set_bands():

    tdb = TDB()

    meg_dat = MegData(tdb, '')

    osc = Osc(default=True)

    meg_dat.set_bands(osc)

def test_cl_set_bands_error():

    tdb = TDB()

    osc = Osc(default=True)

    meg_dat = MegData(tdb, '', osc=osc)

    meg_dat.bands_vertex = True

    with raises(InconsistentDataError):
        meg_dat.set_bands(osc)

def test_cl_import_foof():

    tdb = TDB()

    meg_dat = MegData(tdb, '')

    meg_dat.import_foof('test_v2', get_demo=False, load_type='pickle')

    exp_osc_count = np.array([2, 3])

    assert meg_dat.subnum == 'test_v2'
    assert meg_dat.n_psds == len(meg_dat.slopes) == 2
    assert meg_dat.centers.shape == meg_dat.powers.shape == meg_dat.bws.shape == (2, 8)
    assert np.array_equal(meg_dat.osc_count, exp_osc_count)
    assert meg_dat.has_data

def test_cl_import_foof_error():

    tdb = TDB()

    meg_dat = MegData(tdb, '')
    meg_dat.import_foof('test_v2', get_demo=False, load_type='pickle')

    with raises(InconsistentDataError):
        meg_dat.import_foof('test_v5', get_demo=False, load_type='pickle')

def test_cl_osc_bands_vertex():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.osc_bands_vertex()

    # TODO: ADD MORE PROPER CHECKS HERE!!!
    assert meg_dat.bands_vertex

def tests_cl_osc_bands_vertex_error():

    tdb = TDB()

    meg_dat = MegData(tdb, '')

    with raises(DataNotComputedError):
        meg_dat.osc_bands_vertex()

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

def test_cl_peak_freq_all():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.all_oscs()

    meg_dat.peak_freq(dat='all')

    # TODO: ADD MORE PROPER CHECKS HERE!!! CONSIDER HOW TO USE OSCS FOR PEAKS.
    assert meg_dat.peaks

def test_cl_peak_freq_band():
    pass

def test_cl_peak_freq_errors():

    tdb = TDB()

    meg_dat = MegData(tdb, '')

    with raises(DataNotComputedError):
        meg_dat.peak_freq(dat='all')

    meg_dat.all_osc = True

    with raises(DataNotComputedError):
        meg_dat.peak_freq(dat='all')

def test_cl_calc_osc_param_corrs():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.all_oscs()

    a, b, c = meg_dat.calc_osc_param_corrs()

    # TODO: ADD MORE PROPER CHECKS HERE!!! CONSIDER HOW TO USE OSCS FOR PEAKS.
    assert a.any()

def test_cl_calc_osc_param_corrs_error():

    meg_dat = load_test_meg_subj('test_v2')

    with raises(DataNotComputedError):
        meg_dat.calc_osc_param_corrs()

def test_cl_set_foof_viz():

    meg_dat = load_test_meg_subj('test_v2')

    meg_dat.osc_bands_vertex()

    meg_dat.set_foof_viz()

    assert True
