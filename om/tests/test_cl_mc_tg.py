from __future__ import print_function, division

from types import StringType, ListType
import os

from om.gen import OMDB
import om.cl.mc_tg as mc

############################################################################
################# TESTS - OMEGAMAPPIN - CL_MC_TG - CLASSES #################
############################################################################

def test_mc_tg():
    db = OMDB()

    assert mc.MapCompTG(db)

#############################################################################
############# TESTS - OMEGAMAPPIN - CL_MC_TG - PUBLIC FUNCTIONS #############
#############################################################################

def test_calc_avg_gene_map():
    pass

##############################################################################
############# TESTS - OMEGAMAPPIN - CL_MC_TG - PRIVATE FUNCTIONS #############
##############################################################################

def test_get_map_names_terms():

    db = OMDB()

    names_file = '00-ns_terms.csv'
    names = mc._get_map_names(names_file, db.maps_terms_path)

    assert names
    assert type(names) == ListType
    assert type(names[0]) == StringType

def test_get_map_names_genes():

    db = OMDB()

    names_file = '00-real_gene_names.csv'
    names = mc._get_map_names(names_file, db.maps_genes_path)

    assert names
    assert type(names) is ListType
    assert type(names[0]) is StringType

def test_get_gene_files():

    genes_files_path = mc._get_gene_files('sub1')

    assert genes_files_path
    assert type(genes_files_path) is ListType
    for path in genes_files_path:
        assert os.path.exists(path)

def test_avg_csv_files():
    pass

def test_init_stat_dict():

    assert mc._init_stat_dict(['a', 'b', 'c'])

def test_make_list():
    pass

def test_pull_out_results():
    pass

