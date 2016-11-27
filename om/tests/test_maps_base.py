"""MODULE DOCSTRING - TO FILL IN"""

from om.maps.base import *
from om.maps.base import _init_meg_map_dict
from om.core.db import OMDB
from om.tests.utils import TestDB as TDB

########################################################################################
################# TESTS - OMEGAMAPPIN - CL_MC_BASE - PRIVATE FUNCTIONS #################
########################################################################################

def test_init_meg_map_dict():
    """   """

    bands = ['a', 'b', 'c', 'd']

    assert _init_meg_map_dict(bands)
    assert _init_meg_map_dict(bands, 100)

#######################################################################################
########################## TESTS - OMEGAMAPPIN - MAPS - BASE ##########################
#######################################################################################

def test_mc_base():
    """   """

    db = OMDB()
    assert MapCompBase(db)

def test_cl_load_meg_maps():

    tdb = TDB()
    mc_base = MapCompBase(tdb)

    mc_base.load_meg_maps('test_meg')

    # TODO: ADD FULLER TESTING OF THIS
    assert mc_base.oscs_loaded

def test_cl_load_sl_map():

    tdb = TDB()
    mc_base = MapCompBase(tdb)

    mc_base.load_slope_map('test_slopes')

    # TODO: ADD FULLER TESTING OF THIS
    assert mc_base.slopes_loaded
