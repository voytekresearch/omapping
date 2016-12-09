from __future__ import division

import random
import numpy as np

from sklearn.neighbors import KNeighborsClassifier

from om.core.db import OMDB
from om.meg.single import MegData

##########################################################################################
##########################################################################################
##########################################################################################

def load_subjs(sub_nums, dat_source, db=None):
    """Loads a list of subjects into a list.

    Parameters
    ----------
    sub_nums : list of int
        xx
    dat_source : {'OMEGA', 'HCP'}
        Which database the data comes from.
    db : OMDB() object
        xx

    Returns
    -------
    subjs : list of MegData() objects
        xx
    """

    # Get a database object, unless one is provided
    if not db:
        db = OMDB()

    # Loop through subject numbers to load
    subjs = []
    for sub_num in sub_nums:

        # Load subject of data and process to all oscs
        temp = MegData(db, dat_source)
        temp.import_foof(sub_num)
        temp.all_oscs(verbose=False)

        # Add new subject to output list
        subjs.append(temp)

    return subjs


def test_train_inds(n_oscs, n_train, n_test):
    """

    Parameters
    ----------
    n_oscs : int
        xx
    n_train : list of int
        xx
    n_test : list of int
        xx
    """

    # Check that there are enough oscillations for desired split
    if not n_oscs >= n_train + n_test:
        print('AHHH')

    # Get a shuffled list of possible indices
    inds = range(n_oscs)
    random.shuffle(inds)

    # Pull out random indices for train and test set
    train = inds[0:n_train]
    test = inds[-n_test:]

    return train, test


def check_accuracy_single(dat, lab):
    """Test the accuracy for a list of results, with single answer.

    Parameters
    ----------
    dat : list of int
        xx
    lab : int
        xx

    Returns
    -------
    acc : float
        xx
    """

    #
    n_correct = dat.count(lab)
    acc = n_correct / len(dat)

    return acc

def check_accuracy_all(dat, labels):
    """

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
        print('AHHHHH')

    #
    n_correct = len([True for i, j in zip(dat, labels) if i == j])
    acc = n_correct / len(dat)

    return acc

def knn(subjs, n_train=1500, n_test=50):
    """

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
        tr_inds, te_inds = test_train_inds(subj.n_oscs, n_train, n_test)

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
