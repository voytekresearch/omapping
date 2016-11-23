"""Create fake data for testing OM"""

from __future__ import print_function, division

import os
import sys
import csv
import pickle
import numpy as np

# Import FOOF (use sys to add location to path, then import)
sys.path.append('/Users/thomasdonoghue/Documents/GitCode/omegamappin/')
from om.gen import Osc, save_foof_pickle
from om.tests.helper_test_funcs import TestDB


####

def clear_fake_dat():
    """   """

    tdb = TestDB()

    for key, val in vars(tdb).items():
        if '_path' in key and val:
            _rm_test_files(val)

def _rm_test_files(path):
    """   """

    if os.path.isdir(path):

        files = os.listdir(path)
        for f in files:

            # Skip hidden files
            if f[0] is '.':
                continue

            _rm_test_files(os.path.join(path, f))
    else:

        os.remove(path)

####

def make_fake_foof_dat_1():
    """   """

    tdb = TestDB()

    v1 = (1.0, np.array([5, 10]), np.array([1.0e-22, 1.0e-22]), np.array([1.0, 1.0]))
    v2 = (1.0, np.array([5, 10, 15]), np.array([1.0e-22, 1.0e-22, 1.0e-22]), np.array([1.0, 1.0, 1.0]))

    fake_foof_dat = [v1, v2]

    save_foof_pickle(fake_foof_dat, tdb.foof_path, 'test_v2')

def make_fake_foof_dat_2():
    """   """

    tdb = TestDB()

    v1 = (0.0, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]), np.array([0.5, 1.0, 1.5]))
    v2 = (0.5, np.array([10, 11, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]), np.array([0.5, 1.0, 1.5]))
    v3 = (1.0, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]), np.array([0.5, 1.0, 1.5]))
    v4 = (1.5, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]), np.array([0.5, 1.0, 1.5]))
    v5 = (2.0, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]), np.array([0.5, 1.0, 1.5]))

    fake_foof_dat = [v1, v2, v3, v4, v5]

    save_foof_pickle(fake_foof_dat, tdb.foof_path, 'test_v5')

def make_test_csvs():
    """   """

    tdb = TestDB()

    f_names = ['test1.csv', 'test2.csv']
    f_dat = [[[1.0, 1.0, 1.0, 1.0], [2.0, 2.0, 2.0, 2.0]],
             [[2.0, 2.0, 2.0, 2.0], [3.0, 3.0, 3.0, 3.0]]]

    for ind, f in enumerate(f_names):

        out_file = open(os.path.join(tdb.csvs_path, f), 'wb')
        out_writer = csv.writer(out_file)

        for dat in f_dat[ind]:

            out_writer.writerow(dat)

        out_file.close()

def make_fake_meg_map_data():
    """Create fake MEG map data for testing.

    Test data has MEG aggregate data (prob or score) for 2 bands at 5 vertices.
    """

    tdb = TestDB()

    pickle_name = os.path.join(tdb.maps_path, 'Oscs', 'test_meg.p')

    osc = Osc(input_bands=dict({'a': (3,10), 'b': (10, 40)}))
    test_meg_dat = dict({'a': np.array([1, 1, 1, 1, 1]),
                         'b': np.array([1, 1, 1, 1, 1])})

    dat_out = dict({'bands': osc.bands,
                    'osc_dat': test_meg_dat})

    pickle.dump(dat_out, open(pickle_name, 'wb'))

def make_fake_slope_map_data():
    """Create fake group average slope data for testing.

    Test data has slope data for 5 vertices.
    """

    tdb = TestDB()

    pickle_name = os.path.join(tdb.maps_path, 'Slopes', 'test_slopes.p')

    test_slope_dat = np.array([1, 1, 1, 1, 1])

    dat_out = dict({'slopes': test_slope_dat})

    pickle.dump(dat_out, open(pickle_name, 'wb'))

def make_fake_gene_data():
    """Creates fake genetic data for testing.

    Test data has 5 vertices with 3 genes.
    """

    tdb = TestDB()

    gene_names = ['g-a', 'g-b', 'g-c']

    gene_dat = [[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4], [5, 5, 5]]

    names_f_name = os.path.join(tdb.maps_path, 'Genes', '00-test_gene_names.csv')
    dat_f_name = os.path.join(tdb.maps_path, 'Genes', 'test_gene_estimations',
                              'test_gene_dat_r10.csv')

    _mc_dat(gene_names, gene_dat, names_f_name, dat_f_name)

def make_fake_term_data():
    """Create fake term data for testing.

    Test data has 5 vertices with 2 terms.
    """

    tdb = TestDB()

    term_names = ['t-a', 't-b']

    term_dat = [[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]]

    names_f_name = os.path.join(tdb.maps_path, 'Terms', '00-test_term_names.csv')
    dat_f_name = os.path.join(tdb.maps_path, 'Terms', 'test_term_dat.csv')

    _mc_dat(term_names, term_dat, names_f_name, dat_f_name)

################################################################

def _mc_dat(names, data, name_f_name, dat_f_name):

    names_file = open(name_f_name, 'wb')
    names_writer = csv.writer(names_file)
    names_writer.writerow(names)
    names_file.close()

    dat_file = open(dat_f_name, 'wb')
    dat_writer = csv.writer(dat_file)

    for dat in data:
        dat_writer.writerow(dat)

    dat_file.close()

####################################################################

if __name__ == "__main__":

    clear_fake_dat()
    print('\n\tPrevious test files removed.')

    make_fake_foof_dat_1()
    make_fake_foof_dat_2()

    make_test_csvs()

    make_fake_meg_map_data()
    make_fake_slope_map_data()

    make_fake_gene_data()
    make_fake_term_data()

    print("\tTesting data created.\n")

