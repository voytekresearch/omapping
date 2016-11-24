from __future__ import print_function, division

import os
import csv
import numpy as np
from types import StringType, ListType

from om.gen import OMDB
from helper_test_funcs import TestDB as TDB
import om.cl.mc_tg as mc

####################################################################################
##################### TESTS - OMEGAMAPPIN - CL_MC_TG - CLASSES #####################
####################################################################################

def test_mc_tg():

    db = OMDB()

    assert mc.MapCompTG(db)

###################################################################################
################ TESTS - OMEGAMAPPIN - CL_MC_TG - PUBLIC FUNCTIONS ################
###################################################################################

def test_calc_avg_gene_map():
    pass

####################################################################################
################ TESTS - OMEGAMAPPIN - CL_MC_TG - PRIVATE FUNCTIONS ################
####################################################################################

def test_get_map_names_terms():

    db = OMDB()

    names_file = '00-ns_terms.csv'
    names = mc._get_map_names(names_file, os.path.join(db.maps_path, 'Terms'))

    assert names
    assert type(names) == ListType
    assert type(names[0]) == StringType

def test_get_map_names_genes():

    db = OMDB()

    names_file = '00-real_gene_names.csv'
    names = mc._get_map_names(names_file, os.path.join(db.maps_path, 'Genes'))

    assert names
    assert type(names) is ListType
    assert type(names[0]) is StringType

def test_get_gene_files():

    db = OMDB()

    genes_files_path = mc._get_gene_files('sub1', db)

    assert genes_files_path
    assert type(genes_files_path) is ListType
    for path in genes_files_path:
        assert os.path.exists(path)

def test_avg_csv_files():

    tdb = TDB()

    f_in = [os.path.join(tdb.csvs_path, 'test1.csv'),
            os.path.join(tdb.csvs_path, 'test2.csv')]

    f_out = os.path.join(tdb.csvs_path, 'test_out.csv')

    mc._avg_csv_files(f_in, f_out)

    assert os.path.exists(f_out)

    exp = [[1.5, 1.5, 1.5, 1.5], [2.5, 2.5, 2.5, 2.5]]

    f = open(f_out)
    reader = csv.reader(f)

    for ind, row in enumerate(reader):
        row = [float(i) for i in row]
        assert row == exp[ind]

    f.close()
    os.remove(f_out)

def test_init_stat_dict():

    bands = ['a', 'b', 'c']
    d = mc._init_stat_dict(bands)

    assert d
    bands.append('Slopes')
    assert set(d.keys()) == set(bands)

def test_make_list():
    pass

def test_pull_out_results():
    pass

##############################################################################
############## TESTS - OMEGAMAPPIN - CL_MC_TG - CLASS FUNCTIONS ##############
##############################################################################

def test_load_gene_maps():

    tdb = TDB()

    map_comp = mc.MapCompTG(tdb)

    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')

    # TODO: ADD BETTER TESTING OF THIS
    assert map_comp.genes_loaded
    assert map_comp.gene_maps.shape == (5, 3)


def test_load_term_maps():

    tdb = TDB()

    map_comp = mc.MapCompTG(tdb)

    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')

    # TODO: ADD BETTER TESTING OF THIS
    assert map_comp.terms_loaded
    assert map_comp.term_maps.shape == (5, 2)

def test_calc_corrs_genes_meg_l():

    tdb = TDB()

    map_comp = mc.MapCompTG(tdb)

    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')
    map_comp.load_meg_maps('test_meg')

    for osc in map_comp.bands:
        map_comp.calc_corrs('Genes', osc, method='linear')

    for osc in map_comp.bands:
        assert np.all(map_comp.corrs['Genes'][osc])
        assert np.all(map_comp.p_vals['Genes'][osc])

def test_calc_corrs_genes_slope_l():

    tdb = TDB()

    map_comp = mc.MapCompTG(tdb)

    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')
    map_comp.load_slope_map('test_slopes')

    map_comp.calc_corrs('Genes', 'Slopes', method='linear')

    assert np.all(map_comp.corrs['Genes']['Slopes'])
    assert np.all(map_comp.p_vals['Genes']['Slopes'])

def test_calc_corrs_terms_meg_l():

    tdb = TDB()

    map_comp = mc.MapCompTG(tdb)

    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')
    map_comp.load_meg_maps('test_meg')

    for osc in map_comp.bands:
        map_comp.calc_corrs('Terms', osc, method='linear')

    for osc in map_comp.bands:
        assert np.all(map_comp.corrs['Terms'][osc])
        assert np.all(map_comp.p_vals['Terms'][osc])

def test_calc_corrs_terms_slope_l():

    tdb = TDB()

    map_comp = mc.MapCompTG(tdb)

    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')
    map_comp.load_slope_map('test_slopes')

    map_comp.calc_corrs('Terms', 'Slopes', method='linear')

    assert np.all(map_comp.corrs['Terms']['Slopes'])
    assert np.all(map_comp.p_vals['Terms']['Slopes'])

def test_check_corrs():

    tdb = TDB()

    map_comp = mc.MapCompTG(tdb)

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

def test_unload_data_genes():

    tdb = TDB()

    map_comp = mc.MapCompTG(tdb)

    map_comp.load_gene_maps('test', names_file='00-test_gene_names.csv')
    map_comp.unload_data('Genes')

    assert not map_comp.genes_loaded

def test_unload_data_terms():

    tdb = TDB()

    map_comp = mc.MapCompTG(tdb)

    map_comp.load_term_maps('test_term_dat.csv', names_file='00-test_term_names.csv')
    map_comp.unload_data('Terms')

    assert not map_comp.terms_loaded
