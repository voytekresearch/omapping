# """Tests for OM - meg group."""

# import numpy as np
# from py.test import raises

# from om.meg.group import *
# from om.meg.group import _get_all_osc, _osc_prob, _osc_pow_ratio, _band_sort

# from om.meg.single import MegSubj
# from om.core.db import OMDB
# from om.core.osc import Osc
# from om.tests.utils import TestDB as TDB
# from om.tests.utils import load_test_meg_subj, load_test_meg_gr

# ###################################################################################################
# ###################################################################################################

# def test_group_meg_data():

#     tdb = TDB()
#     osc = Osc()

#     assert MegGroup(tdb, osc)

# def test_add_subject():
#     """

#     TODO: ADD TESTING DIFFERENT ADD SUBJECT SETTINGS
#     """

#     tdb = TDB()
#     osc = Osc(default=True)

#     # Check adding single subject of data
#     meg_group = MegGroup(tdb, osc)
#     meg_subj = load_test_meg_subj('test_v2')
#     meg_group.add_subject(meg_subj)
#     assert meg_group.n_subjs == 1

#     # Test error of adding empty subject of data
#     meg_group = MegGroup(tdb, osc)
#     meg_subj = MegSubj(tdb, '', osc)
#     with raises(DataNotComputedError):
#         meg_group.add_subject(meg_subj)

#     # Test add_vertex_oscs
#     meg_group = MegGroup(tdb, osc)
#     meg_subj = load_test_meg_subj('test_v2')
#     meg_group.add_subject(meg_subj, add_vertex_oscs=True)
#     assert meg_group.has_vertex_oscs

#     # Test add_vertex_exponents
#     meg_group = MegGroup(tdb, osc)
#     meg_subj = load_test_meg_subj('test_v2')
#     meg_group.add_subject(meg_subj, add_vertex_exponents=True)
#     assert meg_group.has_vertex_exponents

#     # Test add_all_oscs
#     meg_group = MegGroup(tdb, osc)
#     meg_subj = load_test_meg_subj('test_v2', all_oscs=True)
#     meg_group.add_subject(meg_subj, add_all_oscs=True)
#     assert meg_group.has_all_osc

#     # Test add_vertex_bands
#     meg_group = MegGroup(tdb, osc)
#     meg_subj = load_test_meg_subj('test_v2', bands_vertex=True)
#     meg_group.add_subject(meg_subj, add_vertex_bands=True)
#     assert meg_group.has_vertex_bands

#     # Test add_peak_freqs
#     meg_group = MegGroup(tdb, osc)
#     meg_subj = load_test_meg_subj('test_v2', all_oscs=True)
#     meg_subj.peak_freq(dat='all')
#     meg_group.add_subject(meg_subj, add_peak_freqs=True)
#     assert meg_group.has_peak_freqs

#     # Test add_demo
#     # TODO


# def test_group_exponent():
#     """   """

#     meg_group = load_test_meg_gr()

#     meg_group.group_exponent()

#     # TODO: ADD MORE TESTING OF THIS
#     assert meg_group.exponent_gr_avg

#     meg_group.group_exponent('median')
#     assert meg_group.exponent_gr_avg

# def test_osc_prob():
#     """   """

#     meg_group = load_test_meg_gr(bands_vertex=True)

#     meg_group.osc_prob()

#     # TODO: ADD MORE TESTING OF THIS
#     assert meg_group.osc_probs
#     assert meg_group.osc_prob_done

# def test_osc_prob_error():
#     """   """

#     meg_group = load_test_meg_gr(bands_vertex=False)

#     with raises(DataNotComputedError):
#         meg_group.osc_prob()

# def test_osc_score():
#     """   """

#     meg_group = load_test_meg_gr(bands_vertex=True)

#     meg_group.osc_prob()
#     meg_group.osc_score()

#     # TODO: ADD MORE TESTING OF THIS
#     assert meg_group.osc_scores
#     assert meg_group.osc_score_done

# def test_osc_map_corrs_prob():
#     """   """

#     meg_group = load_test_meg_gr(bands_vertex=True)

#     meg_group.osc_prob()
#     r, p, b = meg_group.osc_map_corrs('prob')

#     # TODO: ADD TESTING OF THIS
#     assert True

# def test_osc_map_corrs_score():
#     """   """

#     meg_group = load_test_meg_gr(bands_vertex=True, calc_maps=True)

#     r, p, b = meg_group.osc_map_corrs('score')

#     # TODO: ADD TESTING OF THIS
#     assert True

# def test_calc_osc_peak_age():
#     """   """

#     meg_group = load_test_meg_gr(bands_vertex=True, all_osc=True, peaks=True, calc_maps=True)

#     meg_group.age = np.array([20, 20])

#     r, p, b = meg_group.calc_osc_peak_age()

#     # TODO: ADD TESTING OF THIS
#     assert True

# def test_save_gr_exponent():
#     """   """

#     meg_group = load_test_meg_gr()

#     meg_group.group_exponent()

#     meg_group.save_gr_exponent('test_gr_exponent_save')

#     assert True

# def test_save_map():
#     """   """

#     meg_group = load_test_meg_gr(bands_vertex=True, calc_maps=True)

#     meg_group.save_map('prob', 'test_prob_save')
#     meg_group.save_map('score', 'test_score_save')

# def test_set_exponent_viz():
#     """   """

#     meg_group = load_test_meg_gr()

#     meg_group.group_exponent()

#     meg_group.set_exponent_viz()

#     assert True

# def test_set_map_viz():
#     """   """

#     meg_group = load_test_meg_gr(bands_vertex=True, calc_maps=True)

#     meg_group.set_map_viz('prob', 'test_prob_viz_save')
#     meg_group.set_map_viz('score', 'test_score_viz_save')

# #########################################################################################
# ##################### TESTS - OMEGAMAPPIN - MEG - GROUP - FUNCTIONS #####################
# #########################################################################################

# def test_freq_corr_group():
#     """ """

#     meg_group = load_test_meg_gr(vertex_osc=True)

#     f_win = 3
#     corrs, ps, fs = freq_corr_group(meg_group.centers, f_win)

#     assert np.all(fs)

# def test_osc_space_group():
#     """   """
#     pass

# #########################################################################################
# ################# TESTS - OMEGAMAPPIN - MEG - GROUP - PRIVATE FUNCTIONS #################
# #########################################################################################

# def test_get_all_osc():
#     """   """

#     centers = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
#     osc_low = 3
#     osc_high = 7

#     oscs_out = _get_all_osc(centers, osc_low, osc_high)

#     assert len(oscs_out) == 5
#     assert np.all(oscs_out == np.array([3, 4, 5, 6, 7]))

# #def test_osc_prob():
# #    """   """
# #    pass

# def test_osc_pow_ratio():
#     """   """
#     pass

# def test_band_sort():
#     """   """

#     osc = Osc()

#     osc.add_band('b', (12, 14))
#     osc.add_band('a', (4, 5))
#     osc.add_band('c', (15, 19))

#     ord_bands, sort_inds = _band_sort(osc.bands)

#     assert len(ord_bands) == 3
#     assert ord_bands == ['a', 'b', 'c']
#     assert [list(osc.bands.keys())[i] for i in sort_inds] == ord_bands
