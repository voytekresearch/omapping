"""   """

import os
import csv
import numpy as np
import pandas as pd
from types import StringType, ListType
from py.test import raises

from om.maps.tg import *
from om.core.db import OMDB
from om.tests.utils import TestDB as TDB

#####################################################################################
#################### TESTS - OMEGAMAPPIN - MAPS - TG - FUNCTIONS ####################
#####################################################################################

def test_calc_avg_gene_map():
    """   """

    pass
    #tdb = TDB()

    #subj_list = []

    #calc_avg_gene_map(subj_list, 'test_avg')

###################################################################################
######################### TESTS - OMEGAMAPPIN - MAPS - TG #########################
###################################################################################

def test_mc_tg():
    """   """

    db = OMDB()

    assert MapCompTG(db)

def test_load_gene_maps():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')

    # TODO: ADD BETTER TESTING OF THIS
    assert map_comp.genes_loaded
    assert map_comp.gene_maps.shape == (5, 3)

    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')

    assert map_comp.genes_loaded

def test_gene_bad_data():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    with raises(InconsistentDataError):
        map_comp.load_gene_maps('bad_test', names_file='00-test_gene_names.csv')

def test_load_term_maps():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')

    # TODO: ADD BETTER TESTING OF THIS
    assert map_comp.terms_loaded
    assert map_comp.term_maps.shape == (5, 2)

def test_term_bad_data():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    with raises(InconsistentDataError):
        map_comp.load_term_maps('bad_test_term_dat.csv', names_file='00-test_term_names.csv')

def test_calc_corrs_errors():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    with raises(DataNotComputedError):
        map_comp.calc_corrs('Terms', 'a')

    with raises(DataNotComputedError):
        map_comp.calc_corrs('Genes', 'a')

    with raises(UnknownDataTypeError):
        map_comp.calc_corrs('BAD', 'a')

    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')

    with raises(DataNotComputedError):
        map_comp.calc_corrs('Genes', 'Slopes')

    with raises(DataNotComputedError):
        map_comp.calc_corrs('Genes', 'a')

    map_comp.load_meg_maps('test_meg')

    with raises(UnknownDataTypeError):
        map_comp.calc_corrs('Genes', 'BAD')

def test_calc_corrs_genes_meg_l():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')
    map_comp.load_meg_maps('test_meg')

    for osc in map_comp.bands:
        map_comp.calc_corrs('Genes', osc, method='linear')

    for osc in map_comp.bands:
        assert np.all(map_comp.corrs['Genes'][osc])
        assert np.all(map_comp.p_vals['Genes'][osc])

def test_calc_corrs_genes_slope_l():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')
    map_comp.load_slope_map('test_slopes')

    map_comp.calc_corrs('Genes', 'Slopes', method='linear')

    assert np.all(map_comp.corrs['Genes']['Slopes'])
    assert np.all(map_comp.p_vals['Genes']['Slopes'])

def test_calc_corrs_terms_meg_l():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')
    map_comp.load_meg_maps('test_meg')

    for osc in map_comp.bands:
        map_comp.calc_corrs('Terms', osc, method='linear')

    for osc in map_comp.bands:
        assert np.all(map_comp.corrs['Terms'][osc])
        assert np.all(map_comp.p_vals['Terms'][osc])

def test_calc_corrs_terms_slope_l():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')
    map_comp.load_slope_map('test_slopes')

    map_comp.calc_corrs('Terms', 'Slopes', method='linear')

    assert np.all(map_comp.corrs['Terms']['Slopes'])
    assert np.all(map_comp.p_vals['Terms']['Slopes'])

"""
def test_calc_corrs_par():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_meg_maps('test_meg')
    map_comp.load_slope_map('test_slopes')

    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')
    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')

    for osc in map_comp.bands:
        map_comp.calc_corrs('Genes', osc, method='parallel', stop_par=False)

    map_comp.calc_corrs('Terms', 'Slopes', method='parallel', stop_par=True)

    for osc in map_comp.bands:
        assert np.all(map_comp.corrs['Terms'][osc])

    assert np.all(map_comp.corrs['Genes']['Slopes'])
"""

def test_check_corrs():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_meg_maps('test_meg')
    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')
    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')

    for osc in map_comp.bands:
        map_comp.calc_corrs('Genes', osc, method='linear')
        map_comp.calc_corrs('Terms', osc, method='linear')

    for osc in map_comp.bands:
        map_comp.check_corrs('Genes', osc, n_check=2)
        map_comp.check_corrs('Terms', osc, n_check=2)
        map_comp.check_corrs('Genes', osc, n_check=2, top=False)
        map_comp.check_corrs('Terms', osc, n_check=2, top=False)

    assert True

def test_check_corrs_errors():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    with raises(UnknownDataTypeError):
        map_comp.check_corrs('BAD', 'a')

    with raises(DataNotComputedError):
        map_comp.check_corrs('Genes', 'a')

    map_comp.load_meg_maps('test_meg')
    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')

    map_comp.calc_corrs('Terms', 'a')

    with raises(DataNotComputedError):
        map_comp.check_corrs('Terms', 'b')

def test_unload_data_genes():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')
    map_comp.unload_data('Genes')

    assert not map_comp.genes_loaded

def test_unload_data_terms():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')
    map_comp.unload_data('Terms')

    assert not map_comp.terms_loaded

def test_save_corrs():
    """   """

    tdb = TDB()

    map_comp = MapCompTG(tdb)

    map_comp.load_meg_maps('test_meg')
    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')
    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')

    for osc in map_comp.bands:
        map_comp.calc_corrs('Genes', osc, method='linear')
        map_comp.calc_corrs('Terms', osc, method='linear')

    for osc in map_comp.bands:
        map_comp.save_corrs('Genes', osc, 'test')
        map_comp.save_corrs('Terms', osc, 'test')

    with raises(UnknownDataTypeError):
        map_comp.save_corrs('BAD', 'a', 'test')
