from __future__ import print_function, division

import numpy as np
from py.test import raises

#from helper_test_funcs import TestDB as TDB
#from helper_test_funcs import load_test_meg_subj, load_test_meg_gr
from om.tests.utils import TestDB as TDB
from om.tests.utils import load_test_meg_subj, load_test_meg_gr

import om.cl.md_gr as md

#from om.gen import OMDB, Osc
from om.core.db import OMDB
from om.core.osc import Osc

#########################################################################################
###################### TESTS - OMEGAMAPPIN - CL_MD_GROUP - CLASSES ######################
#########################################################################################

def test_group_meg_data():
    """   """

    db = OMDB()
    osc = Osc()

    assert md.GroupMegData(db, osc)

#########################################################################################
##################### TESTS - OMEGAMAPPIN - CL_MD_GROUP - FUNCTIONS #####################
#########################################################################################

def test_get_all_osc():
    """   """

    centers = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    osc_low = 3
    osc_high = 7

    oscs_out = md._get_all_osc(centers, osc_low, osc_high)

    assert len(oscs_out) == 5
    assert np.all(oscs_out == np.array([3, 4, 5, 6, 7]))

def test_osc_prob():
    """   """
    pass


def test_osc_pow_ratio():
    """   """
    pass


def test_band_sort():
    """   """

    osc = Osc()

    osc.add_band('b', [12, 14])
    osc.add_band('a', [4, 5])
    osc.add_band('c', [15, 19])

    ord_bands, sort_inds = md._band_sort(osc.bands)

    assert len(ord_bands) == 3
    assert ord_bands == ['a', 'b', 'c']
    assert [osc.bands.keys()[i] for i in sort_inds]

##########################################################################################
################## TESTS - OMEGAMAPPIN - CL_MD_GROUP - CLASSE FUNCTIONS ##################
##########################################################################################

def test_cl_add_subject():

    db = OMDB()
    osc = Osc()

    meg_group = md.GroupMegData(db, osc)

    meg_subj_dat = load_test_meg_subj('test_v2')

    meg_group.add_subject(meg_subj_dat)

    # TODO: ADD MORE TESTING OF THIS
    assert meg_group

#TODO: ADD TESTING DIFFERENT ADD SUBJECT SETTINGS

def test_cl_group_slope():

    meg_group = load_test_meg_gr()

    meg_group.group_slope()

    # TODO: ADD MORE TESTING OF THIS
    assert meg_group.slopes_gr_avg

    meg_group.group_slope('median')
    assert meg_group.slopes_gr_avg

def test_cl_osc_prob():

    meg_group = load_test_meg_gr(bands_vertex=True)

    meg_group.osc_prob()

    # TODO: ADD MORE TESTING OF THIS
    assert meg_group.osc_probs
    assert meg_group.osc_prob_done

def test_osc_prob_error():

    meg_group = load_test_meg_gr(bands_vertex=False)

    with raises(md.DataNotComputedError):
        meg_group.osc_prob()

def test_cl_osc_score():

    meg_group = load_test_meg_gr(bands_vertex=True)

    meg_group.osc_prob()
    meg_group.osc_score()

    # TODO: ADD MORE TESTING OF THIS
    assert meg_group.osc_scores
    assert meg_group.osc_score_done

def test_cl_osc_map_corrs_prob():

    meg_group = load_test_meg_gr(bands_vertex=True)

    meg_group.osc_prob()
    r, p, b = meg_group.osc_map_corrs('prob')

    # TODO: ADD TESTING OF THIS
    assert True

def test_cl_osc_map_corrs_score():

    meg_group = load_test_meg_gr(bands_vertex=True, calc_maps=True)

    r, p, b = meg_group.osc_map_corrs('score')

    # TODO: ADD TESTING OF THIS
    assert True

def test_cl_calc_osc_peak_age():

    meg_group = load_test_meg_gr(bands_vertex=True, all_osc=True, peaks=True, calc_maps=True)

    meg_group.age = np.array([20, 20])

    r, p, b = meg_group.calc_osc_peak_age()

    # TODO: ADD TESTING OF THIS
    assert True

def test_cl_freq_corr():

    meg_group = load_test_meg_gr(vertex_osc=True)

    r, p = meg_group.freq_corr(1)

    # TODO: ADD TESTING OF THIS
    assert True

def test_save_gr_slope():

    meg_group = load_test_meg_gr()

    meg_group.group_slope()

    meg_group.save_gr_slope('test_gr_slope_save')

    assert True

def test_save_map():

    meg_group = load_test_meg_gr(bands_vertex=True, calc_maps=True)

    meg_group.save_map('prob', 'test_prob_save')
    meg_group.save_map('score', 'test_score_save')

def test_set_slope_viz():

    meg_group = load_test_meg_gr()

    meg_group.group_slope()

    meg_group.set_slope_viz()

    assert True

def test_set_map_viz():

    meg_group = load_test_meg_gr(bands_vertex=True, calc_maps=True)

    meg_group.set_map_viz('prob', 'test_prob_viz_save')
    meg_group.set_map_viz('score', 'test_score_viz_save')
