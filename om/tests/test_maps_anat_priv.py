"""   """

import numpy as np

import om.maps.anat as anat

######################################################################################
################### TESTS - OMEGAMAPPIN - ANAT - PRIVATE FUNCTIONS ###################
######################################################################################

def test_extract_lr_anat():
    """   """

    anat_labels = ['left_caudal_anterior_cingulate',
                   'right_caudal_middle_frontal']

    labels, lrs = anat._extract_lr(anat_labels, 'anat')

    assert len(labels) == 2
    assert len(lrs) == 2
    assert labels[0] == 'caudalanteriorcingulate'
    assert labels[1] == 'caudalmiddlefrontal'
    assert lrs[0] is 'L'
    assert lrs[1] is 'R'

def test_extract_lr_elec():
    """   """

    elec_labels = ['parsopercularis R', 'postcentral L',]

    labels, lrs = anat._extract_lr(elec_labels, 'elec')

    assert len(labels) == 2
    assert len(lrs) == 2
    assert labels[0] == 'parsopercularis'
    assert labels[1] == 'postcentral'
    assert lrs[0] is 'R'
    assert lrs[1] is 'L'

def test_clean_label_anat():
    """   """

    label_1 = 'left_caudal_anterior_cingulate'
    label_2 = 'right_caudal_anterior_cingulate'
    lr_1 = 'L'
    lr_2 = 'R'

    label_1 = anat._clean_label(label_1, lr_1, 'anat')
    label_2 = anat._clean_label(label_2, lr_2, 'anat')

    assert label_1 == 'caudalanteriorcingulate'
    assert label_2 == 'caudalanteriorcingulate'

def test_clean_label_elec():
    """   """

    label = 'parsopercularis R'
    lr = 'R'

    label = anat._clean_label(label, lr, 'elec')

    assert label == 'parsopercularis'

def test_mat_mult():
    """   """

    dat = np.array([1, 2, 3])

    res = anat._mat_mult(dat)
    exp = np.array([[1, 2, 3], [2, 4, 6], [3, 6, 9]])

    assert res.shape == (3, 3)
    assert np.all(res == exp)
