"""MODULE DOCSTRING"""
from __future__ import print_function

import os
import csv
import numpy as np
from scipy.stats.stats import pearsonr

from om.core.io import load_meg_list
from om.core.db import OMDB, check_db
from om.core.osc import Osc
from om.meg.single import MegData

##########################################################################################
##########################################################################################
##########################################################################################

def get_twin_data(db=None, file_name='00-HCP_Subjects_RESTRICTED.csv'):
    """Extract twin status data from data file.

    Parameters
    ----------
    db : OMDB() object, optional
        xx
    f_name : str, optional
        xx

    Returns
    -------
    mz_twins : 2d array
        IDs for MZ twins and their parents - columns: [twin_id, mother_id, father_id].
    dz_twins : 2d array
        Ds for DZ twins and their parents - columns: [twin_id, mother_id, father_id].
    twin_list : list of int
        Subject IDs for all twins.
    not_twin_list : list of int
        Subject IDs for non-twin subjects.

    Notes
    -----
    TODO:
    """

    # Check db, initialize if not provided
    db = check_db(db)

    # Set up file indices
    id_ind = 0
    twin_ind = 2
    zyg_ind = 3
    moth_ind = 4
    fath_ind = 5

    # Get and sort twins, by type
    twin_list = []
    mz_twins = np.empty(shape=[0, 3], dtype='int_')
    dz_twins = np.empty(shape=[0, 3], dtype='int_')
    not_twin_list = []

    # Open file, and use csv reader for parsing
    with open(os.path.join(db.meg_path, file_name)) as f_name:
        reader = csv.reader(f_name, delimiter=',')

        # Loop through each row in the file
        for row in reader:

            # If subject is a twin, add to running list of all twins
            if row[twin_ind] == 'Twin':

                twin_list.append(int(row[id_ind]))

                # If MZ, add to running list of MZ twins
                if row[zyg_ind] == 'MZ':

                    mz_twins = np.vstack((mz_twins, [int(row[id_ind]), int(row[moth_ind]),
                                                     int(row[fath_ind])]))

                # If NotMZ, add to running list of DZ twins
                elif row[zyg_ind] == 'NotMZ':

                    dz_twins = np.vstack((dz_twins, [int(row[id_ind]), int(row[moth_ind]),
                                                     int(row[fath_ind])]))

            # If not a twin, add to running list of not twins
            elif row[twin_ind] == 'NotTwin':
                not_twin_list.append(int(row[id_ind]))

    return mz_twins, dz_twins, twin_list, not_twin_list


def match_twins(dat, parent_ind=1):
    """Match twin pairs.

    Parameters
    ----------
    dat : 2d array
        IDs for twins and their parents - columns: [twin_id, mother_id, father_id].
    parent_ind : {1, 2}, optional (default=1)
        Which parent to use to match twins, which should be irrelevant. Defaults to mother.

    Returns
    -------
    twin_pairs : list of tuple of (int, int)
        Each element within the list contains the subject IDs for a twin pair.
    single_twins : list of tuple of (int)
        Each element within the list contains the ID of a twins that is unmatched.

    Notes
    -----
    TODO: EXPLAIN LOGIC.
    """

    # Pull out relevant data from input matrix
    all_parent_ids = dat[:, parent_ind]
    all_subj_ids = dat[:, 0]

    # Get all the unique parent IDs
    unique_parents = set(list(all_parent_ids))

    # Initliaze variables to store data
    twin_pairs = []
    single_twins = []

    # For a given parent, find their kids
    for parent in unique_parents:

        # Find all kids for current parents
        kids = list(np.where(all_parent_ids == parent)[0])

        # If only one subject found, add to unpaired twins
        if len(kids) == 1:
            single_twins.append(tuple(all_subj_ids[kids]))

        # If two subjects found, add as a pair of twins
        elif len(kids) == 2:
            twin_pairs.append(tuple(all_subj_ids[kids]))
        else:
            print('AHHHHH')

    return twin_pairs, single_twins


def check_complete_pairs(twin_ids, available_files):
    """Check which twin pairs have both sets of subject data available.

    Parameters
    ----------
    twin_ids : list of tuple of (int, int)
        xx
    available_files : list of int
        xx

    Returns
    -------
    complete_pairs : ?
        xx
    """

    complete_pairs = []

    for pair in twin_ids:

        if pair[0] in available_files and pair[1] in available_files:
            complete_pairs.append(pair)

    return complete_pairs


def rm_twin_pairs(all_pairs, twin_pairs):
    """Given all possible subject pairs, remove twins leaving only unrelated pairs.

    Parameters
    ----------
    all_pairs : list of tuple of (int, int)
        xx
    twin_pairs : list of tuple of (int, int)
        xx

    Returns
    -------
    non_twins : list of tuple of (int, int)
        xx
    """

    non_twins = list(set(all_pairs) - set(twin_pairs))

    return non_twins


def compare_pair(pair_inds, osc=None, db=None, dat_source='HCP'):
    """Compares center frequency data for a pairing of MEG subjects.

    Parameters
    ----------
    pair_inds : list of int
        xx
    osc : Osc() object, optional
        xx
    db : OMDB() object, optional
        xx
    dat_source : {'HCP', 'OMEGA'}, optional (default='HCP')
        Which database to use for data.

    Returns
    -------
    corr_dat : 2d array, shape=(n_bands, 2), row = [r-val, p-val]
        Results of the correlation within each band between subjects.
    """

    # Check db, initialize if not provided
    db = check_db(db)

    # Load data
    #dat = load_pair(pair_inds, osc_bands_vert=True, osc=osc, db=db, dat_source=dat_source)
    dat = load_meg_list(pair_inds, osc_bands_vert=True, osc=osc, db=db, dat_source=dat_source)


    # Initialize to store correlation results
    corr_dat = np.zeros([4, 2])

    # Compare center frequencies within oscillatory bands
    for ind, band in enumerate(osc.bands):
        corr_dat[ind, 0], corr_dat[ind, 1] = pearsonr(dat[0].oscs[band][:, 0],
                                                      dat[1].oscs[band][:, 0])

    return corr_dat


def compare_slope(pair_inds, db=None, dat_source='HCP'):
    """Compares slope value data for a pair of MEG subjects.

    Parameters
    ----------
    pair_inds : ?
        xx
    db : OMDB() object, optional
        xx
    dat_source : {'HCP', 'OMEGA'}, optional (default='HCP')
        Which database to use for data.

    Returns
    -------
    corr_dat : 1d array, [r-val, p-val]
        Results of the correlation between subjects.
    """

    # Check db, initialize if not provided
    db = check_db(db)

    # Load data
    #dat = load_pair(pair_inds, db=db, dat_source=dat_source)
    dat = load_meg_list(pair_inds, db=db, dat_source=dat_source)

    # Initialize to store correlation results
    corr_dat = np.zeros(2)

    # Compare slopes
    corr_dat[0], corr_dat[1] = pearsonr(dat[0].slopes, dat[1].slopes)

    return corr_dat


def print_twin_results(corr_dat, labels):
    """Print out correlation results.

    Parameters
    ----------
    corr_dat : ?
        xx
    labels : ?
        xx
    """

    for ind, label in enumerate(labels):
        print('\t', label, '\t : ', '{:5.4f}'.format(corr_dat[ind, 0]))
