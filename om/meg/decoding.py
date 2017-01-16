"""Decoding subject ID from magnetophysiological data."""

from __future__ import division

import random
import numpy as np

from sklearn.neighbors import KNeighborsClassifier

from om.core.db import OMDB, check_db
from om.meg.single import MegSubj
from om.core.errors import InconsistentDataError

##########################################################################################
##########################################################################################
##########################################################################################

def arrange_dat(subjs):
    """Arrange data from subjects to be able to run decoding on them.

    Parameters
    ----------
    subjs : ?
        xx

    Returns
    -------
    all_dat : 2d array
        xx
    all_labels : 1d array
        xx
    """

    # Initilialize data variables
    all_dat = np.empty((0, 3))
    all_labels = np.array([], dtype=int)

    # Loop through each subject
    for ind, subj in enumerate(subjs):

        # Pull out data from current subject
        all_dat = np.vstack([all_dat, np.array([subj.centers_all, subj.powers_all, subj.bws_all]).T])

        # Set the train labels for current subject
        all_labels = np.append(all_labels, np.ones(len(subj.centers_all), dtype=int) * ind)

    return all_dat, all_labels


def split_inds(n_oscs, n_train, n_test):
    """NOTE: OLD! Use 'train_test_split' from sklearn

    Split indices into train and test groups.

    Parameters
    ----------
    n_oscs : int
        The number of oscillations that can be split.
    n_train : list of int
        Number of data points for training group.
    n_test : list of int
        Number of data points for testing group.

    Returns
    -------
    train : list of int
        Indices to be used for training data.
    test : list of int
        Indices to be used for test data.
    """

    # Check that there are enough oscillations for desired split
    if not n_oscs >= n_train + n_test:
        raise InconsistentDataError('Number of data points is insufficient for requested split.')

    # Get a shuffled list of possible indices
    inds = range(n_oscs)
    random.shuffle(inds)

    # Pull out random indices for train and test set
    train = inds[0:n_train]
    test = inds[-n_test:]

    return train, test


def check_accuracy_single(dat, label):
    """NOTE: OLD! UNUSED. Use sklearn funcs.

    Test the accuracy for a list of results, with single answer.

    Parameters
    ----------
    dat : list of int
        Predicted labels of a set of test data.
    label : int
        The correct label.

    Returns
    -------
    acc : float
        Prediction accuracy for given data.
    """

    # Check how accuracy of all results
    n_correct = dat.count(label)
    acc = n_correct / len(dat)

    return acc


def check_accuracy_all(dat, labels):
    """NOTE: OLD! UNUSED. Use sklearn funcs.

    Parameters
    ----------
    dat : list of int
        xx
    labels : list of int
        xx

    Returns
    -------
    acc : float
        xx
    """

    # Check data and labels are the same length
    if not len(dat) == len(labels):
        raise InconsistentDataError('Length of data does not match labels.')

    #
    n_correct = len([True for i, j in zip(dat, labels) if i == j])
    acc = n_correct / len(dat)

    return acc


def knn(subjs, n_train=1500, n_test=50):
    """NOTE: OLD! UNUSED. Use sklearn funcs.

    Parameters
    ----------
    subjs : ?
        xx
    n_train : int
        xx
    n_test : int
        xx
    """

    # Initilialize data variables
    all_dat = np.empty((0, 3))
    all_labels = np.array([], dtype=int)
    all_test_inds = []
    all_test_labels = []

    # Loop through each subject
    for ind, subj in enumerate(subjs):

        # Get random sets of indices to pull out data
        tr_inds, te_inds = split_inds(subj.n_oscs, n_train, n_test)

        # Pull out the training data from current subject
        all_dat = np.vstack([all_dat, np.array([subj.centers_all[tr_inds], subj.powers_all[tr_inds],
                                                subj.bws_all[tr_inds]]).T])

        # Set the train labels for current subject
        all_labels = np.append(all_labels, np.ones(n_train, dtype=int) * ind)

        # Pull out the test indices and test labels for current subject
        all_test_inds.append(te_inds)
        temp_labels = [ind] * n_test
        all_test_labels = all_test_labels + temp_labels

    # Iniliaize KNN model and fit training data
    neigh = KNeighborsClassifier(n_neighbors=5)
    neigh.fit(all_dat, all_labels)

    # Test model
    results = []
    for s_ind, t_inds in enumerate(all_test_inds):

        for t_i in t_inds:

            test_dat = np.array([subjs[s_ind].centers_all[t_i], subjs[s_ind].powers_all[t_i],
                                 subjs[s_ind].bws_all[t_i]], ndmin=2)
            results.append(neigh.predict(test_dat)[0])

    acc = check_accuracy_all(results, all_test_labels)

    return acc
