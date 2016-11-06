from __future__ import print_function, division
from om.gen import OMDB
import om.cl.mc_anat as mc

######################################################################
################## TESTS - OMEGAMAPPIN - CL_MC_ANAT ##################
######################################################################


def test_mc_anat():
    """   """

    db = OMDB()

    assert mc.MapCompAnat(db)


def test_mat_mult():
    """   """

    pass


def test_make_list():
    """   """

    pass


def test_pull_out_results():
    """   """

    pass