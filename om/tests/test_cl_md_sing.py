from __future__ import print_function, division
import numpy as np
import om.cl.md as md
from om.gen import Osc

#################################################################################
########################## TESTS - OMEGAMAPPIN - CL_MD ##########################
#################################################################################

def test_get_osc():
    """   """
    pass


def test_get_all_osc():
    """   """
    pass


def test_get_single_osc_power():
    """   """
    pass

def test_get_demo_csv():
    """   """
    pass


def test_osc_prob():
    """   """
    pass


def test_osc_pow_ratio():
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


def test_band_sort():
    """   """

    # Initialize osc object, add some bands
    #osc_bands = Osc()
    #osc_bands.add_band('b', [12, 14])
    #osc_bands.add_band('a')
    pass




