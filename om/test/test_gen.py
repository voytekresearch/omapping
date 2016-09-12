from om.gen import *

###############################################
########## TESTS - OMEGAMAPPIN - GEN ##########
###############################################

def test_clean_file_list():
    """Test that clean_file_list works properly. 

    Tests that the function returns one item to chosen query, 
        and that it is case-insensitive. 
    """

    # Initiate vars to pass in
    files = ['not_this_one.py', 'or_this.txt', 'FiNd_WorD.yes']
    string = 'find_word'

    # Run clean file list
    out = clean_file_list(files, string)

    # Check that returns a list of 1 item, that is 'FiNd_WorD.ys'
    assert(len(out) == 1)
    assert(out[0] is files[2])


def test_extract_psd():
    """

    TODO: Finish this test. 


    # Initialize input vars
    psd = []
    freqs = range(16)
    f_low = 10
    f_high = 5

    # Run extract psd
    psd_out, freqs_out = extract_psd(psd, freqs, f_low, f_high)

    # Check if answers as expected
    assert()
    """
    pass

def test_get_sub_nums_first():
    """   """

    subj_files = ['123_test.py', '234_test.py', '345_test.py']
    out = get_sub_nums(subj_files, 'first')

    assert(len(out) == len(subj_files))
    assert(out[0] == 123)
    assert(out[-1] == 345)


def test_get_sub_nums_last():
    """   """

    subj_files = ['test_123.py', 'test_234.py', 'test_345.py']
    out = get_sub_nums(subj_files, 'last')

    assert(len(out) == len(subj_files))
    assert(out[0] == 123)
    assert(out[-1] == 345)


def test_get_cur_subj():
    """   """

    cur_subj = 123
    subj_files = ['123_test.npz', '456_meg.csv', '12_3.aa', 'aa_231.ok']
    out = get_cur_subj(cur_subj, subj_files)

    assert(str(cur_subj) in out)

