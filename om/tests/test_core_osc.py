"""   """

from py.test import raises

from om.core.osc import *

#########################################################################################
########################## TESTS - OMEGAMAPPIN - CLASSES - OSC ##########################
#########################################################################################

def test_osc():
    """Test that Osc() returns properly, with different inputs."""

    assert Osc()
    assert Osc(default=True)
    assert Osc(input_bands=dict({'Test Band': (5, 10)}))

def test_add_band():
    """Test that Osc.add_band() adds band properly."""

    osc = Osc()

    osc.add_band('test', [0, 100])

    assert 'test' in osc.bands.keys()
    assert osc.bands['test'][0] == 0
    assert osc.bands['test'][1] == 100

def test_osc_add_band_error():
    """Test that Osc.add_band() returns errors appropriately."""

    osc = Osc()

    with raises(InconsistentDataError):
        osc.add_band('bad', [0, 0])

    with raises(InconsistentDataError):
        osc.add_band('worse', [10, 5])

def test_rm_band():
    """Test that Osc.rm_band() removes bands properly."""

    osc = Osc()

    osc.add_band('test', [1, 2])

    osc.rm_band('test')

    with raises(KeyError):
        assert osc.bands['test']
