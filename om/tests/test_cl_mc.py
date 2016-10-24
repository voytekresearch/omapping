from __future__ import print_function, division
import om.cl.mc as mc

#####################################################################
#################### TESTS - OMEGAMAPPIN - CL_MC ####################
#####################################################################

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
