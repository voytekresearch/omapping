"""   """

import os
from py.test import raises

from om.core.utils import *
from om.tests.utils import TestDB as TDB

##########################################################################################
##########################################################################################
##########################################################################################

def test_clean_file_list():
    """Test that clean_file_list() works properly.

    Tests that the function returns one item to chosen query,
        and that it is case-insensitive.
    """

    files = ['not_this_one.py', 'or_this.txt', 'FiNd_WorD.yes']
    string = 'find_word'

    out = clean_file_list(files, string)

    assert len(out) == 1
    assert out[0] is files[2]

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

    n_rois = 6
    roi_lr = ['L', 'L', 'L', 'R', 'R', 'R']

    st_x, en_x, st_y, en_y = get_section('all', n_rois, roi_lr)
    assert st_x == st_y == 0
    assert en_x == en_y == n_rois

    st_x, en_x, st_y, en_y = get_section('left', n_rois, roi_lr)
    assert st_x == st_y == 0
    assert en_x == en_y == 2

    st_x, en_x, st_y, en_y = get_section('right', n_rois, roi_lr)
    assert st_x == st_y == 3
    assert en_x == en_y == 5

    st_x, en_x, st_y, en_y = get_section('lr', n_rois, roi_lr)
    assert st_x == 0; assert st_y == 3
    assert en_x == 2; assert en_y == 5

    st_x, en_x, st_y, en_y = get_section('rl', n_rois, roi_lr)
    assert st_x == 3; assert st_y == 0
    assert en_x == 5; assert en_y == 2

    with raises(InconsistentDataError):
        get_section('bad', n_rois, roi_lr)

def test_extract_psd():
    """Test that extract_psd() works properly."""

    psd = np.zeros([5, 16])
    freqs = np.array(range(16))
    f_low = 5
    f_high = 10

    psd_out, freqs_out = extract_psd(psd, freqs, f_low, f_high)
    n_row, n_col = np.shape(psd_out)

    assert freqs_out.min() == f_low
    assert freqs_out.max() == f_high
    assert (n_row == 5) & (n_col == len(freqs_out))

def test_avg_csv_files():
    """   """

    tdb = TDB()

    f_in = [os.path.join(tdb.csvs_path, 'test1.csv'),
            os.path.join(tdb.csvs_path, 'test2.csv')]

    f_out = os.path.join(tdb.csvs_path, 'test_out.csv')

    avg_csv_files(f_in, f_out)

    assert os.path.exists(f_out)

    exp = [[1.5, 1.5, 1.5, 1.5], [2.5, 2.5, 2.5, 2.5]]

    f = open(f_out)
    reader = csv.reader(f)

    for ind, row in enumerate(reader):
        row = [float(i) for i in row]
        assert row == exp[ind]

    f.close()
    os.remove(f_out)

    avg_csv_files(f_in, f_out, 'median')
    assert os.path.exists(f_out)
    os.remove(f_out)
