from __future__ import print_function, division
import numpy as np
import om.cl.md_sing as md
from om.gen import OMDB

##################################################################################
################### TESTS - OMEGAMAPPIN - CL_MD_SING - CLASSES ###################
##################################################################################

def test_meg_data():
    """   """

    db = OMDB()

    assert md.MegData(db)


##########################################################################################
################## TESTS - OMEGAMAPPIN - CL_MD_SING - PRIVATE FUNCTIONS ##################
##########################################################################################

def test_get_single_osc():
    """   """
    pass

def test_get_single_osc_power():
    """   """
    pass

def test_get_demo_csv():
    """   """
    pass

def test_osc_peak():
    """   """

    # Initialize data
    centers = np.array([1, 2, 3, 4])
    osc_low = 1.5
    osc_high = 3.5

    # Run function
    test_ans = md._osc_peak(centers, osc_low, osc_high)

    # Check the returned answer is close to the actual answer
    actual_ans = (2 + 3) / 2
    assert np.isclose(test_ans, actual_ans)

##########################################################################################
################### TESTS - OMEGAMAPPIN - CL_MD_SING - CLASS FUNCTIONS ###################
##########################################################################################

def test_cl_osc_bands_vertex():
    pass

def test_cl_all_oscs():
    pass

def test_cl_peak_freq():
    pass

def test_cl_calc_osc_param_corrs():
    pass
