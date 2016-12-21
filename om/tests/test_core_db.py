"""Testing of database functions and classes for OM."""

import os
from types import ListType

from om.core.db import OMDB

#########################################################################################
########################## TESTS - OMEGAMAPPIN - CLASSES - OMDB #########################
#########################################################################################

def test_omdb():
    """Test that OMDB() returns properly."""

    assert OMDB(auto_gen=False)

def test_omdb_gen_paths():
    """Test that gen_paths method of OMDB works."""

    db = OMDB(auto_gen=False)
    db.gen_paths()

    assert db

def test_omdb_paths():
    """Test that all defined OMDB paths exist."""

    db = OMDB()
    db.gen_paths()

    # Loops through all paths, checking they exist
    #  Skips vars with '_path' marker, and empty variables
    for key, val in vars(db).items():
        if '_path' in key and val:
            assert os.path.exists(val)

def test_check_dat_files_psd():
    """Test that OMDB() check_dat_files method runs properly, for PSD files."""

    db = OMDB()

    sub_nums, source = db.check_dat_files('PSD', verbose=True)

    assert sub_nums
    assert source
    assert len(sub_nums) == len(source)

    sub_nums, source = db.check_dat_files('PSD', 'OMEGA', verbose=True)
    sub_nums, source = db.check_dat_files('PSD', 'HCP', verbose=True)

def test_check_dat_files_foof():
    """Test that OMDB() check_dat_files method runs properly, for FOOF files."""

    db = OMDB()

    sub_nums, source = db.check_dat_files('foof', verbose=True)

    assert type(sub_nums) == type(source) == ListType
    assert len(sub_nums) == len(source)

    sub_nums, source = db.check_dat_files('foof', 'OMEGA', verbose=True)
    sub_nums, source = db.check_dat_files('foof', 'HCP', verbose=True)

def test_check_res_files_meg():
    """Test that OMDB() check_res_files method runs properly, for MEG files."""

    db = OMDB()

    files = db.check_res_files('meg', verbose=True)

    assert type(files) == ListType

def test_check_res_files_maps():
    """Test that OMDB() check_res_files method runs properly, for MAPS files."""

    db = OMDB()

    files = db.check_res_files('maps', verbose=True)

    assert type(files) == ListType

def test_check_map_files():
    """Test that OMDB() check_maps_files method runs properly."""

    db = OMDB()

    osc, sl, gene, term = db.check_map_files(verbose=True, return_files=True)

    assert type(osc) == type(sl) == type(gene) == type(term) == ListType
