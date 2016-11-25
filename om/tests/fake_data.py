"""Create fake data for testing OM"""

from __future__ import print_function, division

import os
import sys
import csv
import pickle
import numpy as np
import scipy.io as sio

# Import FOOF (use sys to add location to path, then import)
sys.path.append('/Users/thomasdonoghue/Documents/GitCode/omegamappin/')
from om.gen import Osc, save_foof_pickle
from om.tests.helper_test_funcs import TestDB

####################
####################
####################

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

def make_test_file_directory(base_path):
    """   """

    # Corrs Data
    os.mkdir('Corrs')
    cor_data = ['Genes', 'Terms']
    cor_data_type = ['csv', 'npz']

    for cor_dat in cor_data:
        os.mkdir('Corrs/' + cor_dat)

        for cor_dat_type in cor_data_type:
            os.mkdir('Corrs' + cor_dat + cor_dat_type)

    # Maps Data
    os.mkdir('Maps')
    maps_data = ['Genes', 'Oscs', 'Slopes', 'Terms']
    for maps_dat in maps_data:
        os.mkdir('Maps/' + maps_dat)

    # MEG Data
    os.mkdir('MEG')

    #
    os.mkdir('csvs')
    os.mkdir('foof')

################################
################################
################################

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

    v1 = (0.0, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]),
          np.array([0.5, 1.0, 1.5]))
    v2 = (0.5, np.array([10, 11, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]),
          np.array([0.5, 1.0, 1.5]))
    v3 = (1.0, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]),
          np.array([0.5, 1.0, 1.5]))
    v4 = (1.5, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]),
          np.array([0.5, 1.0, 1.5]))
    v5 = (2.0, np.array([4, 12, 30]), np.array([1.0e-22, 1.5e-22, 2.0e-22]),
          np.array([0.5, 1.0, 1.5]))

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

    _mc_dat(gene_dat, dat_f_name, gene_names, names_f_name)

def make_bad_gene_data():

    tdb = TestDB()

    bad_gene_dat = [[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]]

    bad_dat_f_name = os.path.join(tdb.maps_path, 'Genes', 'bad_test_gene_estimations',
                              'bad_test_gene_dat_r10.csv')

    _mc_dat(bad_gene_dat, bad_dat_f_name)

def make_fake_term_data():
    """Create fake term data for testing.

    Test data has 5 vertices with 2 terms.
    """

    tdb = TestDB()

    term_names = ['t-a', 't-b']

    term_dat = [[1, 1], [2, 2], [3, 3], [4, 4], [5, 5]]

    names_f_name = os.path.join(tdb.maps_path, 'Terms', '00-test_term_names.csv')
    dat_f_name = os.path.join(tdb.maps_path, 'Terms', 'test_term_dat.csv')

    _mc_dat(term_dat, dat_f_name, term_names, names_f_name)

def make_bad_term_data():

    tdb = TestDB()

    bad_term_dat = [[1], [2], [3], [4], [5]]

    bad_dat_f_name = os.path.join(tdb.maps_path, 'Terms', 'bad_test_term_dat.csv')

    _mc_dat(bad_term_dat, bad_dat_f_name)

def make_fake_anat_data():

    tdb = TestDB()

    roi_labels = np.array([['left_test'], ['right_test']], dtype=np.object)
    connectivity = np.array([[0.5, 0.5], [0.5, 0.5]])

    sio.savemat(os.path.join(tdb.maps_path, 'Anat', 'test_anat.mat'),
                dict({'roi_labels': roi_labels,
                      'connectivity': connectivity}))

def make_fake_scout_data():

    tdb = TestDB()

    dt = [('Vertices', 'O'), ('Seed', 'f4'), ('Color', 'O'), ('Label', 'S10'),
          ('Function', 'S10'), ('Region', 'S10'), ('Handles', 'S10')]
    scout_dat = np.zeros((2,), dtype=dt)

    scout_dat[0]['Vertices'] = [1, 2, 5]
    scout_dat[0]['Seed'] = 1111
    scout_dat[0]['Color'] = [1, 1, 1]
    scout_dat[0]['Label'] = 'test L'
    scout_dat[0]['Function'] = 1111
    scout_dat[0]['Region'] = 'Mean'
    scout_dat[0]['Handles'] = ''

    scout_dat[1]['Vertices'] = [3, 4]
    scout_dat[1]['Seed'] = 1111
    scout_dat[1]['Color'] = [1, 1, 1]
    scout_dat[1]['Label'] = 'test R'
    scout_dat[1]['Function'] = 1111
    scout_dat[1]['Region'] = 'Mean'
    scout_dat[1]['Handles'] = ''

    sio.savemat(os.path.join(tdb.maps_path, 'Scouts', 'test_scout.mat'),
                {'Scouts': scout_dat})

##############################################################
##############################################################
##############################################################

def _mc_dat(data, dat_f_name, names=None, name_f_name=None):

    if names:
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

    make_bad_gene_data()
    make_bad_term_data()

    make_fake_anat_data()
    make_fake_scout_data()

    print("\tTesting data created.\n")

