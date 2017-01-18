"""Heritability analysis, using the twin data from the HCP data."""
from __future__ import print_function, division

import os
import csv
import numpy as np
from scipy.stats.stats import pearsonr

from om.core.io import load_meg_list
from om.core.db import OMDB, check_db
from om.core.osc import Osc, check_bands, CheckBands
from om.meg.single import MegSubj

#TODO: Make check_bands functionality into a decorator?

##################################################################################################
##################################################################################################
##################################################################################################

def get_twin_data(db=None, file_name='00-HCP_Subjects_RESTRICTED.csv'):
    """Extract twin status data from data file.

    Parameters
    ----------
    db : OMDB() object, optional
        Database object for omegamappin project.
    f_name : str, optional
        File name for the metadata file to use.

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
        List of subject IDs for twin pairs.
    single_twins : list of tuple of (int)
        IDs for twins who's pair is not available in MEG data.

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
            # TODO: fix this up
            print('AHHHHH')

    return twin_pairs, single_twins


def check_complete_pairs(twin_ids, available_files):
    """Check which twin pairs have both sets of subject data available.

    Parameters
    ----------
    twin_ids : list of tuple of (int, int)
        Twin ID pairs.
    available_files : list of int
        All subject data files that are available.

    Returns
    -------
    complete_pairs : list of tuple of (int, int)
        Twin ID pairs for which both subject's data are available.
    """

    # Loop through all twin pairs
    complete_pairs = []
    for pair in twin_ids:

        # Keep pair if both IDs are available
        if pair[0] in available_files and pair[1] in available_files:
            complete_pairs.append(pair)

    return complete_pairs


def rm_twin_pairs(all_pairs, twin_pairs):
    """Given all possible subject pairs, remove twins leaving only unrelated pairs.

    Parameters
    ----------
    all_pairs : list of tuple of (int, int)
        All possible subject ID pairs.
    twin_pairs : list of tuple of (int, int)
        Twin ID pairs.

    Returns
    -------
    list of tuple of (int, int)
        Subject ID pairs for all non-twins.
    """

    # Remove all pairs that are twin pairs
    return list(set(all_pairs) - set(twin_pairs))


def comp_peak_freq(pairs_dat, peak_type='band', avg='median'):
    """Compares the peak frequencies, within oscillatory bands, between pairs of subjects.

    Parameters
    ----------
    pairs_dat : list of list of [MegSubj, MegSubj]
        List of subject pairs to compare.
    peak_type : {'bands', 'all'}
        What data to use to compute the peak frequency.
    avg : {'median', 'mean'}
        Which type of averaging to use.

    Returns
    -------
    peak_avg : 1d array
        Average peak frequency differences, length = n_bands.
    peak_dat : 2d array
        All data for all subjects, [n_subjects, n_bands].
    """

    # Compare peak frequency between pairs and take average across all pairs
    peak_dat = np.asarray([_comp_pf_pair(pair, peak_type=peak_type, avg=avg) for pair in pairs_dat])
    peak_avg = np.mean(peak_dat, axis=0)

    return peak_avg, peak_dat


def comp_osc_space(pairs_dat):
    """Compares the spatial overlap of oscillation bands between pairs of subjects.

    Parameters
    ----------
    pairs_dat : list of list of [MegSubj, MegSubj]
        List of subject pairs to compare.

    Returns
    -------
    space_avg : 1d array
        Average percent spatial overlap across all pairs.
    space_dat : 2d array
        Spatial overlap results for all pairs.
    """

    # Compare spatial overlap of oscillation bands between pairs and take average across all pairs
    space_dat = np.asarray([_comp_space_pair(pair) for pair in pairs_dat])
    space_avg = np.mean(space_dat, axis=0)

    return space_avg, space_dat


def comp_osc_param(pairs_dat, osc_param):
    """Compares oscillatory band parameters between pairs of subjects.

    Parameters
    ----------
    pairs_dat : list of list of [MegSubj, MegSubj]
        List of subject pairs to compare.

    Returns
    -------
    param_avg : 1d array
        Average correlation between oscillation parameters across all pairs.
    param_dat : 2d array
        Oscillation parameter correlation data for all pairs.
    """

    # Compare oscillatory band parameters between pairs and take average across all pairs
    param_dat = np.asarray([_comp_osc_pair(pair, osc_param) for pair in pairs_dat])
    param_avg = np.mean(param_dat, axis=0)

    return param_avg, param_dat


def comp_slope(pairs_dat):
    """Compares slope values between pairs of subjects.

    Parameters
    ----------
    pairs_dat : list of list of [MegSubj, MegSubj]
        List of subject pairs to compare.

    Returns
    -------
    slope_avg : 1d array
        Average slope corrlation across all pairs.
    slope_dat : 2d array
        Slope correlation data data for all pairs.
    """

    # Compare slope values between pairs and take average across all pairs
    slope_dat = np.asarray([_comp_sl_pair(pair) for pair in pairs_dat])
    slope_avg = np.mean(slope_dat, axis=0, keepdims=True)

    return slope_avg, slope_dat


def print_twin_results_corr(corr_dat, labels):
    """Print out correlation results.

    Parameters
    ----------
    corr_dat : 2d array
        Matrix of correlation data to print out.
    labels : list of str
        Labels for what data each row corresponds to.
    """

    for ind, label in enumerate(labels):
        print('\t', label, '\t : ', '{:5.4f}'.format(corr_dat[ind, 0]))


def print_twin_results_vec(vec_dat, labels):
    """Print out comparison results, stored in 1d vector.

    Parameters
    ----------
    vec_dat : 1d array
        Vector of data to print out.
    labels : list of str
        Labels for what data each row corresponds to.
    """

    for ind, label in enumerate(labels):
        print('\t', label, '\t : ', '{:5.4f}'.format(vec_dat[ind]))

###################################################################################################
###################################################################################################
###################################################################################################

@CheckBands
def _comp_pf_pair(dat, bands, peak_type='band', avg='median'):
    """Compare peak frequencies between a pair of MEG subjects.

    Parameters
    ----------
    dat : list of [MegSubj, MegSubj]
        Pair of subjects to compare.
    bands : dict()
        Oscillation band definitions.
    peak_type : {'bands', 'all'}, optional
        What data to use to compute the peak frequency.
    avg : {'median', 'mean'}, optional
        Which type of averaging to use.

    Returns
    -------
    peak_dists : 1d array
        Distance (in Hz) between oscillatory band peak frequencies of pair.
    """

    # Calculate peak freqs for subjects
    for subj in dat:
        subj.peak_freq(peak_type, avg)

    peak_dists = np.zeros([len(bands)])

    # Calculate distance between peak freqs
    for ind, band in enumerate(bands):
        peak_dists[ind] = abs(dat[1].peaks[band] - dat[0].peaks[band])

    return peak_dists


@CheckBands
def _comp_space_pair(dat, bands):
    """Compare the spatial overlap of oscillatory bands between two subjects.

    Parameters
    ----------
    dat : list of [MegSubj, MegSubj]
        Pair of subjects to compare.
    bands : dict()
        Oscillation band definitions.

    Returns
    -------
    res : 1d array
        Vector of spatial overlap, in percent, for each oscillatory band.

    NOTE:
    - Update hard coded value for number of vertices.
        - Add to MegSubj object as an attribute?
    """

    # Initialize variable to store output data
    res = np.zeros([len(bands)])

    # Loop through each oscillatory band
    for ind, band in enumerate(bands):

        # Create boolean arrays of vertices with oscillation, for each subject, then compare
        #   Use the number of vertices that are the same (have / don't have osc) divided by n_verts
        res[ind] = sum((dat[0].oscs[band][:, 0] > 0) == (dat[1].oscs[band][:, 0] > 0)) / 7500

    return res


@CheckBands
def _comp_osc_pair(dat, osc_param, bands):
    """Compares oscillation band parameters for a pair of MEG subjects.

    Parameters
    ----------
    dat : list of [MegSubj, MegSubj]
        Pair of subjects to compare.
    osc_param : {0, 1, 2}
        Which oscillatory parameter to compare: {0: CF, 1: Power, 2: BW}.
    bands : dict()
        Oscillation band definitions.

    Returns
    -------
    corr_dat : 2d array, shape=(n_bands, 2), row = [r-val, p-val]
        Results of the correlation within each band between subjects.
    """

    # Initialize to store correlation results
    corr_dat = np.zeros([len(bands), 2])

    # Compare center frequencies within oscillatory bands
    for ind, band in enumerate(bands):
        corr_dat[ind, 0], corr_dat[ind, 1] = pearsonr(dat[0].oscs[band][:, osc_param],
                                                      dat[1].oscs[band][:, osc_param])

    return corr_dat


def _comp_sl_pair(dat):
    """Compares slope value data for a pair of MEG subjects.

    Parameters
    ----------
    dat : list of [MegSubj, MegSubj]
        Pair of subjects to compare.

    Returns
    -------
    corr_dat : 1d array, [r-val, p-val]
        Results of the correlation between subjects.
    """

    # Initialize to store correlation results
    corr_dat = np.zeros(2)

    # Compare slopes
    corr_dat[0], corr_dat[1] = pearsonr(dat[0].slopes, dat[1].slopes)

    return corr_dat
