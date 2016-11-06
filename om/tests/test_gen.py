"""MODULE DOCSTRING - TO FILL IN
"""

import numpy as np
from py.test import raises
from om.gen import *

###############################################################################
########################## TESTS - OMEGAMAPPIN - OMDB #########################
###############################################################################

def test_omdb():
    """Test that OMDB() returns properly."""

    # Check that OMDB returns properly, with different data source settings
    assert OMDB()
    assert OMDB('OMEGA')
    assert OMDB('HCP')
    with raises(UnknownDataSourceError):
        assert OMDB('bad data selection')

###############################################################################
########################## TESTS - OMEGAMAPPIN - OSC ##########################
###############################################################################

def test_osc():
    """Test that Osc() returns properly, with different inputs."""

    # Check tha Osc returns properly
    assert Osc()
    assert Osc(default=True)
    assert Osc(input_bands=dict({'Test Band': (5, 10)}))


def test_add_band():
    """Test that Osc.add_band() adds band properly."""

    # Initialize an Osc object
    osc = Osc()

    # Add a new band to osc
    osc.add_band('test', [0, 100])

    # Check it added properly
    assert 'test' in osc.bands.keys()
    assert osc.bands['test'][0] == 0
    assert osc.bands['test'][1] == 100


def test_osc_add_band_error():
    """Test that Osc.add_band() returns errors appropriately."""

    # Initialize an Osc object
    osc = Osc()

    # Test that it raises an error given illegal band definitions
    with raises(InconsistentDataError):
        osc.add_band('bad', [0, 0])

    with raises(InconsistentDataError):
        osc.add_band('worse', [10, 5])


def test_rm_band():
    """Test that Osc.rm_band() removes bands properly."""

    # Initialize an Osc object
    osc = Osc()

    # Add a band
    osc.add_band('test', [1, 2])

    # Remove the band
    osc.rm_band('test')

    # Check that band is no longer there
    with raises(KeyError):
        assert osc.bands['test']

###############################################################################
########################## TESTS - OMEGAMAPPIN - OSC ##########################
###############################################################################

def test_fig_info():
    """Test that FigInfo() returns properly."""

    # Check tha Osc returns properly
    assert FigInfo()

###################################################################################
############################ TESTS - OMEGAMAPPIN - GEN ############################
###################################################################################

def test_clean_file_list():
    """Test that clean_file_list() works properly.

    Tests that the function returns one item to chosen query,
        and that it is case-insensitive.
    """

    # Initiate vars to pass in
    files = ['not_this_one.py', 'or_this.txt', 'FiNd_WorD.yes']
    string = 'find_word'

    # Run clean file list
    out = clean_file_list(files, string)

    # Check that returns a list of 1 item, that is 'FiNd_WorD.ys'
    assert len(out) == 1
    assert out[0] is files[2]


def test_extract_psd():
    """Test that extract_psd() works properly."""

    # Initialize input vars
    psd = np.zeros([5, 16])
    freqs = np.array(range(16))
    f_low = 5
    f_high = 10

    # Run extract psd
    psd_out, freqs_out = extract_psd(psd, freqs, f_low, f_high)
    n_row, n_col = np.shape(psd_out)

    # Check if answers as expected
    assert freqs_out.min() == f_low
    assert freqs_out.max() == f_high
    assert (n_row == 5) & (n_col == len(freqs_out))


def test_get_sub_nums_first():
    """Test that get_sub_num() works properly with numbers first."""

    subj_files = ['123_test.py', '234_test.py', '345_test.py']
    out = get_sub_nums(subj_files, 'first')

    assert len(out) == len(subj_files)
    assert out[0] == 123
    assert out[-1] == 345


def test_get_sub_nums_last():
    """Test that get_sub_num() works properly with numbers last."""

    subj_files = ['test_123.py', 'test_234.py', 'test_345.py']
    out = get_sub_nums(subj_files, 'last')

    assert len(out) == len(subj_files)
    assert out[0] == 123
    assert out[-1] == 345


def test_get_cur_subj():
    """Tests that get_cur_subj() returns correct file name."""

    cur_subj = 123
    subj_files = ['123_test.npz', '456_meg.csv', '12_3.aa', 'aa_231.ok']
    out = get_cur_subj(cur_subj, subj_files)

    assert str(cur_subj) in out


def test_rm_files_ext():
    """   """

    subj_files = ['123_test.npz', '456_meg.csv', '12_3.aa', 'aa_231.ok']
    files_out = rm_files_ext(subj_files)

    assert len(files_out) == 4
    assert '123_test' in files_out[0]
    assert not 'npz' in files_out[0]
    assert not '.ok' in files_out[3]


def test_get_section():
    """   """
    pass
