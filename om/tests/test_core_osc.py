"""Tests for OM - osc class & functions."""

from py.test import raises

from om.core.osc import *
from om.core.osc import _check_band

###################################################################################################
###################################################################################################

def test_osc():
    """Test that Osc() returns properly, with different inputs."""

    assert Osc()
    assert Osc(default=True)
    assert Osc(input_bands=dict({'Test Band': (5, 10)}))

def test_add_band():
    """Test that Osc.add_band() adds band properly."""

    osc = Osc()

    osc.add_band('test', (0, 100))

    assert 'test' in osc.bands.keys()
    assert osc.bands['test'][0] == 0
    assert osc.bands['test'][1] == 100

def test_rm_band():
    """Test that Osc.rm_band() removes bands properly."""

    osc = Osc()

    osc.add_band('test', (1, 2))

    osc.rm_band('test')

    with raises(KeyError):
        assert osc.bands['test']

###################################################################################################
###################################################################################################

def test_check_band():
    """Test that _check_band returns errors appropriately."""

    # Check for error checking that key is a string
    with raises(InconsistentDataError):
        _check_band(111, (5, 10))

    # Check for error checking that limits have length 2
    with raises(InconsistentDataError):
        _check_band('string', (5, 10, 15))

    # Check for proper band definitions & order - 1
    with raises(InconsistentDataError):
        _check_band('bad', (0, 0))

    # Check for proper band definitions & order - 2
    with raises(InconsistentDataError):
        _check_band('worse', (10, 5))
