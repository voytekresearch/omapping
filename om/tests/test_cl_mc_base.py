"""MODULE DOCSTRING - TO FILL IN
"""

from __future__ import print_function, division

from om.gen import OMDB
import om.cl.mc_base as mc

########################################################################################
###################### TESTS - OMEGAMAPPIN - CL_MC_BASE - CLASSES ######################
########################################################################################

def test_mc_base():
    """   """

    db = OMDB()
    assert mc.MapCompBase(db)

########################################################################################
################# TESTS - OMEGAMAPPIN - CL_MC_BASE - PRIVATE FUNCTIONS #################
########################################################################################

def test_init_meg_map_dict():
    """   """

    bands = ['a', 'b', 'c', 'd']

    assert mc._init_meg_map_dict(bands)
    assert mc._init_meg_map_dict(bands, 100)

########################################################################################
################## TESTS - OMEGAMAPPIN - CL_MC_BASE - CLASS FUNCTIONS ##################
########################################################################################

def test_cl_load_meg_maps():
    pass

def test_cl_load_sl_map():
    pass
