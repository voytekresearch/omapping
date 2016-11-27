"""   """

import numpy as np

import om.meg.group as group
from om.core.osc import Osc

#########################################################################################
################# TESTS - OMEGAMAPPIN - MEG - GROUP - PRIVATE FUNCTIONS #################
#########################################################################################

def test_get_all_osc():
    """   """

    centers = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    osc_low = 3
    osc_high = 7

    oscs_out = group._get_all_osc(centers, osc_low, osc_high)

    assert len(oscs_out) == 5
    assert np.all(oscs_out == np.array([3, 4, 5, 6, 7]))

def test_osc_prob():
    """   """
    pass


def test_osc_pow_ratio():
    """   """
    pass


def test_band_sort():
    """   """

    osc = Osc()

    osc.add_band('b', [12, 14])
    osc.add_band('a', [4, 5])
    osc.add_band('c', [15, 19])

    ord_bands, sort_inds = group._band_sort(osc.bands)

    assert len(ord_bands) == 3
    assert ord_bands == ['a', 'b', 'c']
    assert [osc.bands.keys()[i] for i in sort_inds]
