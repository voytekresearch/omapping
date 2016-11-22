"""MODULE DOCSTRING - TO FILL IN
"""

import numpy as np
import os
from types import ListType
from py.test import raises

import om.gen as gen

#########################################################################################
########################## TESTS - OMEGAMAPPIN - CLASSES - OMDB #########################
#########################################################################################

def test_omdb():
    """Test that OMDB() returns properly."""

    assert gen.OMDB()
    assert gen.OMDB('OMEGA')
    assert gen.OMDB('HCP')
    with raises(gen.UnknownDataSourceError):
        assert gen.OMDB('bad data selection')

def test_omdb_paths():
    """Test that all defined OMDB paths exist."""

    db = gen.OMDB()

    # Loops through all paths, checking they exist
    #  Skips vars with '_path' marker, and empty variables
    for key, val in vars(db).items():
        if '_path' in key and val:
            assert os.path.exists(val)

def test_check_dat_files_psd():

    db = gen.OMDB()

    sub_nums, source = db.check_dat_files('PSD', verbose=False)

    assert sub_nums
    assert source
    assert len(sub_nums) == len(source)

def test_check_dat_files_foof():

    db = gen.OMDB()

    sub_nums, source = db.check_dat_files('foof', verbose=False)

    assert type(sub_nums) == type(source) == ListType
    assert len(sub_nums) == len(source)

def test_check_res_files_md():

    db = gen.OMDB()

    files = db.check_res_files('md', verbose=False)

    assert type(files) == ListType

def test_check_res_files_mc():

    db = gen.OMDB()

    files = db.check_res_files('mc', verbose=False)

    assert type(files) == ListType

def test_check_map_files():

    db = gen.OMDB()

    osc, sl, gene, term = db.check_map_files(verbose=False, return_files=True)

    assert type(osc) == type(sl) == type(gene) == type(term) == ListType

#########################################################################################
########################## TESTS - OMEGAMAPPIN - CLASSES - OSC ##########################
#########################################################################################

def test_osc():
    """Test that Osc() returns properly, with different inputs."""

    assert gen.Osc()
    assert gen.Osc(default=True)
    assert gen.Osc(input_bands=dict({'Test Band': (5, 10)}))

def test_add_band():
    """Test that Osc.add_band() adds band properly."""

    osc = gen.Osc()

    osc.add_band('test', [0, 100])

    assert 'test' in osc.bands.keys()
    assert osc.bands['test'][0] == 0
    assert osc.bands['test'][1] == 100

def test_osc_add_band_error():
    """Test that Osc.add_band() returns errors appropriately."""

    osc = gen.Osc()

    with raises(gen.InconsistentDataError):
        osc.add_band('bad', [0, 0])

    with raises(gen.InconsistentDataError):
        osc.add_band('worse', [10, 5])

def test_rm_band():
    """Test that Osc.rm_band() removes bands properly."""

    osc = gen.Osc()

    osc.add_band('test', [1, 2])

    osc.rm_band('test')

    with raises(KeyError):
        assert osc.bands['test']

#########################################################################################
######################### TESTS - OMEGAMAPPIN - CLASSES - FIGINFO ########################
#########################################################################################

def test_fig_info():
    """Test that FigInfo() returns properly."""

    assert gen.FigInfo()

##########################################################################################
###################### TESTS - OMEGAMAPPIN - GEN - PUBLIC FUNCTIONS ######################
##########################################################################################

def test_clean_file_list():
    """Test that clean_file_list() works properly.

    Tests that the function returns one item to chosen query,
        and that it is case-insensitive.
    """

    files = ['not_this_one.py', 'or_this.txt', 'FiNd_WorD.yes']
    string = 'find_word'

    out = gen.clean_file_list(files, string)

    assert len(out) == 1
    assert out[0] is files[2]

def test_extract_psd():
    """Test that extract_psd() works properly."""

    psd = np.zeros([5, 16])
    freqs = np.array(range(16))
    f_low = 5
    f_high = 10

    psd_out, freqs_out = gen.extract_psd(psd, freqs, f_low, f_high)
    n_row, n_col = np.shape(psd_out)

    assert freqs_out.min() == f_low
    assert freqs_out.max() == f_high
    assert (n_row == 5) & (n_col == len(freqs_out))

def test_get_sub_nums_first():
    """Test that get_sub_num() works properly with numbers first."""

    subj_files = ['123_test.py', '234_test.py', '345_test.py']
    out = gen.get_sub_nums(subj_files, 'first')

    assert len(out) == len(subj_files)
    assert out[0] == 123
    assert out[-1] == 345

def test_get_sub_nums_last():
    """Test that get_sub_num() works properly with numbers last."""

    subj_files = ['test_123.py', 'test_234.py', 'test_345.py']
    out = gen.get_sub_nums(subj_files, 'last')

    assert len(out) == len(subj_files)
    assert out[0] == 123
    assert out[-1] == 345

def test_get_cur_subj():
    """Tests that get_cur_subj() returns correct file name."""

    cur_subj = 123
    subj_files = ['123_test.npz', '456_meg.csv', '12_3.aa', 'aa_231.ok']
    out = gen.get_cur_subj(cur_subj, subj_files)

    assert str(cur_subj) in out

def test_rm_files_ext():
    """   """

    subj_files = ['123_test.npz', '456_meg.csv', '12_3.aa', 'aa_231.ok']
    files_out = gen.rm_files_ext(subj_files)

    assert len(files_out) == 4
    assert '123_test' in files_out[0]
    assert not 'npz' in files_out[0]
    assert not '.ok' in files_out[3]


def test_get_section():
    """   """
    pass

###########################################################################################
###################### TESTS - OMEGAMAPPIN - GEN - PRIVATE FUNCTIONS ######################
###########################################################################################

def test_check_files():
    pass

def test_find_last():

    l = [1, 2, 3, 1, 2]
    w = 1

    out = gen._find_last(l, w)

    assert out == 3
