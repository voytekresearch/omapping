from __future__ import print_function, division
import numpy as np
import om.cl.md_gr as md
from om.gen import OMDB, Osc

#################################################################################
################## TESTS - OMEGAMAPPIN - CL_MD_GROUP - CLASSES ##################
#################################################################################

def test_group_meg_data():
    """   """

    db = OMDB()
    osc = Osc()

    assert md.GroupMegData(db, osc)

###################################################################################
################## TESTS - OMEGAMAPPIN - CL_MD_GROUP - FUNCTIONS ##################
###################################################################################

def test_get_all_osc():
    """   """
    pass


def test_osc_prob():
    """   """
    pass


def test_osc_pow_ratio():
    """   """
    pass


def test_band_sort():
    """   """

    # Initialize osc object, add some bands
    #osc_bands = Osc()
    #osc_bands.add_band('b', [12, 14])
    #osc_bands.add_band('a')
    pass
