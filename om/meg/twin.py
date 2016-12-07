"""MODULE DOCSTRING"""

import os
import csv
import numpy as np

from om.core.db import OMDB

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
