from __future__ import print_function, division
from om.gen import OMDB
import om.cl.mc_tg as mc

######################################################################
################### TESTS - OMEGAMAPPIN - CL_MC_TG ###################
######################################################################


def test_mc_tg():
    """   """

    db = OMDB()

    assert mc.MapCompTG(db)


def test_get_map_names():
    """   """

    pass


def test_init_stat_dict():
    """   """

    assert mc._init_stat_dict(['a', 'b', 'c'])