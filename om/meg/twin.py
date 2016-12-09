"""MODULE DOCSTRING"""

import os
import csv
import numpy as np
from scipy.stats.stats import pearsonr

from om.core.db import OMDB
from om.core.osc import Osc
from om.meg.single import MegData

##########################################################################################
##########################################################################################
##########################################################################################

def get_twin_data():
    """Extract twin status data from data file.

    Returns
    -------
    mz_twins : ?
        xx
    dz_twins : ?
        xx
    twin_list : ?
        xx
    not_twin_list : ?
        xx
    """

    # Initialize database object, and set file name
    db = OMDB()
    file_name = '00-HCP_Subjects_RESTRICTED.csv'

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


def match_twins(dat, parent_ind):
    """Match twin pairs.

    Parameters
    ----------
    dat : ?
        xx
    parent_ind : ?
        xx

    Returns
    -------
    twin_pairs : ?
        xx
    single_twins : ?
        xx
    """

    # Pull out relevant data from input matrix
    all_parents = dat[:, parent_ind]
    all_ids = dat[:, 0]

    #
    unique_parents = set(list(all_parents))

    # Initliaze variables to store data
    pair_inds = []
    twin_pairs = []
    single_inds = []
    single_twins = []

    #
    for parent in unique_parents:

        check_pair = list(np.where(all_parents == parent)[0])

        #
        if len(check_pair) == 1:

            single_inds.append(check_pair)
            single_twins.append(list(all_ids[check_pair]))

        #
        elif len(check_pair) == 2:

            pair_inds.append(check_pair)
            twin_pairs.append(list(all_ids[check_pair]))

    return twin_pairs, single_twins


def check_complete_pairs(twin_ids, available_files):
    """Check which twin pairs have both sets of subject data available.

    Parameters
    ----------
    twin_ids : ?
        xx
    available_files : ?
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
    all_pairs : list of int
        xx
    twin_pairs : list of int
        xx

    Returns
    -------
    all_pairs : list of int
        xx
    """

    for pair in all_pairs:

        for twins in twin_pairs:

            if set(pair) == set(twins):
                all_pairs.remove(pair)

    return all_pairs


def compare_pair(pair_inds, db=None):
    """Compares center frequency data for a pairing of MEG subjects.

    Parameters
    ----------
    pair_inds : list of int
        xx

    Returns
    -------
    corr_dat : 2d array
        xx
    """

    # Initialize database object, unless one is supplied
    if not db:
        db = OMDB()

    # Set up oscillation band definition, and dat source
    osc = Osc(default=True)
    dat_source = 'HCP'

    # Initialize data object and load data - pair data-A
    pair_a = MegData(db, dat_source, osc)
    pair_a.import_foof(pair_inds[0], get_demo=False)
    pair_a.osc_bands_vertex()

    # Initialize data object and load data - pair data-B
    pair_b = MegData(db, dat_source, osc)
    pair_b.import_foof(pair_inds[1], get_demo=False)
    pair_b.osc_bands_vertex()

    # Initialize to store correlation results
    corr_dat = np.zeros([4, 2])

    # Compare center frequencies within oscillatory bands
    for ind, band in enumerate(osc.bands):
        corr_dat[ind, 0], corr_dat[ind, 1] = pearsonr(pair_a.oscs[band][:, 0],
                                                      pair_b.oscs[band][:, 0])

    return corr_dat
