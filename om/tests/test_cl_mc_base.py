from __future__ import print_function, division
from om.gen import OMDB
import om.cl.mc_base as mc

######################################################################
################## TESTS - OMEGAMAPPIN - CL_MC_BASE ##################
######################################################################

def test_mc_base():
    """   """

    db = OMDB()
    assert MapCompBase(db)

def test_get_map_names():
    """   """

    pass

def test_init_meg_map_dict():
    """   """

    bands = ['a', 'b', 'c', 'd']

    assert mc._init_meg_map_dict(bands)
    assert mc._init_meg_map_dict(bands, 100)

def test_init_stat_dict():
    """   """

    assert mc._init_stat_dict(['a', 'b', 'c'])

def test_mat_mult():
    """   """

    pass

def test_make_list():
    """   """

    pass

def test_pull_out_results():
    """   """

    pass
