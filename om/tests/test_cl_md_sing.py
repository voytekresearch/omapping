from __future__ import print_function, division

import numpy as np

from helper_test_funcs import load_test_meg_subj
from helper_test_funcs import TestDB as TDB
import om.cl.md_sing as md
from om.gen import OMDB, Osc

##################################################################################
################### TESTS - OMEGAMAPPIN - CL_MD_SING - CLASSES ###################
##################################################################################

def test_meg_data():
    """   """

    db = OMDB()

    assert md.MegData(db, '')


##########################################################################################
################## TESTS - OMEGAMAPPIN - CL_MD_SING - PRIVATE FUNCTIONS ##################
##########################################################################################

def test_get_single_osc():
    """   """
    pass

def test_get_single_osc_power():
    """   """
    pass

def test_get_demo_csv():
    """   """
    pass

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

def test_cl_import_foof():

    tdb = TDB()

    meg_dat = md.MegData(tdb, '')

    meg_dat.import_foof('test1', get_demo=False, load_type='pickle')

    exp_osc_count = np.array([2, 3])

    assert meg_dat.subnum == 'test1'
    assert meg_dat.n_psds == len(meg_dat.slopes) == 2
    assert meg_dat.centers.shape == meg_dat.powers.shape == meg_dat.bws.shape == (2, 8)
    assert np.array_equal(meg_dat.osc_count, exp_osc_count)
    assert meg_dat.has_data

def test_cl_osc_bands_vertex():

    meg_dat = load_test_meg_subj('test1')

    meg_dat.osc_bands_vertex()

    # TODO: ADD MORE PROPER CHECKS HERE!!!
    assert meg_dat.bands_vertex

def test_cl_all_oscs():

    meg_dat = load_test_meg_subj('test1')

    meg_dat.all_oscs(verbose=False)

    # TODO: ADD MORE PROPER CHECKS HERE!!!
    assert meg_dat.all_osc

def test_cl_peak_freq():

    meg_dat = load_test_meg_subj('test1')

    meg_dat.all_oscs()

    meg_dat.peak_freq()

    # TODO: ADD MORE PROPER CHECKS HERE!!! CONSIDER HOW TO USE OSCS FOR PEAKS.
    assert meg_dat.peaks

def test_cl_calc_osc_param_corrs():

    meg_dat = load_test_meg_subj('test1')

    meg_dat.all_oscs()

    a, b, c = meg_dat.calc_osc_param_corrs()

    # TODO: ADD MORE PROPER CHECKS HERE!!! CONSIDER HOW TO USE OSCS FOR PEAKS.
    assert a.any()
