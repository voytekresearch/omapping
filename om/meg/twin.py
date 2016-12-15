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
    twin_pairs : list of list of int
        Each list within the list contains the subject IDs for a twin pair.
    single_twins : list of list of int
        Each list within the list contains the ID of a twins that is unmatched.

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
    #pair_inds = []
    #single_inds = []

    # For a given parent, find their kids
    for parent in unique_parents:

        # Find all kids for current parents
        kids = list(np.where(all_parent_ids == parent)[0])

        # If only one subject found, add to unpaired twins
        if len(kids) == 1:
            single_twins.append(list(all_subj_ids[kids]))
            #single_inds.append(kids)

        # If two subjects found, add as a pair of twins
        elif len(kids) == 2:
            twin_pairs.append(list(all_subj_ids[kids]))
            #pair_inds.append(kids)
        else:
            print('AHHHHH')

    return twin_pairs, single_twins


def check_complete_pairs(twin_ids, available_files):
    """Check which twin pairs have both sets of subject data available.

    Parameters
    ----------
    twin_ids : ?
        xx
    available_files : ?
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


def load_pair(pair_inds, osc_bands_vert=False, osc=None, db=None, dat_source='HCP'):
    """Load a pair of MEG subjects into MegData objects.

    Parameters
    ----------
    pair_inds : ?
        xx
    osc_band_vert : boolean, optional (default:False)
        xx
    osc : Osc() object, optional
        xx
    db : OMDB() object, optional
        xx
    dat_source : {'HCP', 'OMEGA'}, optional (default='HCP')
        Which database to use for data.

    Returns
    -------
    dat_out : list of MegData() objects
        xx
    """

    # Initialize database object, unless one is supplied
    if not db:
        db = OMDB()

    # Initialize data object to return
    dat_out = []

    # Loop through subject pair, loading and processing data
    for subj_id in pair_inds:

        # Initialize MegData, load foof dat
        temp = MegData(db, dat_source, osc)
        temp.import_foof(subj_id, get_demo=False)

        # Convert to oscillation vertex bands, if requested
        if osc_bands_vert:
            temp.osc_bands_vertex()

        # Add subject to output object
        dat_out.append(temp)

    return dat_out


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
    corr_dat : 2d array
        xx
    """

    # Initialize database object, unless one is supplied
    if not db:
        db = OMDB()

    # Load data
    dat = load_pair(pair_inds, osc_bands_vert=True, osc=osc, db=db, dat_source=dat_source)

    # Initialize to store correlation results
    corr_dat = np.zeros([4, 2])

    # Compare center frequencies within oscillatory bands
    for ind, band in enumerate(osc.bands):
        corr_dat[ind, 0], corr_dat[ind, 1] = pearsonr(dat[0].oscs[band][:, 0],
                                                      dat[1].oscs[band][:, 0])

    return corr_dat


def compare_pair_old(pair_inds, db=None):
    """Compares center frequency data for a pair of MEG subjects.

    Parameters
    ----------
    pair_inds : list of int
        xx
    db : OMDB() object, optional
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
    corr_dat : ?
        xx
    """

    # Initialize database object, unless one is supplied
    if not db:
        db = OMDB()

    # Load data
    dat = load_pair(pair_inds, db=db, dat_source=dat_source)

    # Initialize to store correlation results
    corr_dat = np.zeros(2)

    # Compare slopes
    corr_dat[0], corr_dat[1] = pearsonr(dat[0].slopes, dat[1].slopes)

    return corr_dat
