"""Load and save functions for the OM project."""

import os
import csv
import pickle
import datetime

import scipy.io as sio

from om.meg.single import MegSubj
from om.core.db import check_db
from om.core.utils import clean_file_list, get_cur_subj
from om.core.errors import UnknownDataSourceError, UnknownDataTypeError

###################################################################################################
###################################################################################################

def load_meg_psds(dat_source, meg_path, subj_num):
    """Loads a requested subject's PSD-MEG data.

    Parameters
    ----------
    dat_source : {'OMEGA', 'HCP'}
        Which database subject comes from.
    meg_path : str
        Path where data to load is located.
    subj_num : int
        Subject identifier number.

    Returns
    -------
    psd : 2d array
        Matrix of PSD results for each vertex [n_verts, n_freqs].
    freqs : 1d array
        Vector of the frequencies estimated in the PSD.
    """

    # Set file name and get full path
    if dat_source is 'OMEGA':
        mat_file = 'psd_source_median_' + str(subj_num)
    elif dat_source is 'HCP':
        mat_file = 'PSD_Source_' + str(subj_num)
    file_name = os.path.join(meg_path, dat_source, ('Subject_' + str(subj_num)), mat_file)

    # Load MEG PSD data from matfile
    data_mat = sio.loadmat(file_name, appendmat=True, struct_as_record=False, squeeze_me=True)

    # Pull out data from dictionary
    freqs = data_mat['Freqs']
    psd = data_mat['TF']

    return psd, freqs


def save_obj_pickle(obj, dat_type, save_name, db=None):
    """Save a custom object out to a pickle file.

    Parameters
    ----------
    obj : obj()
        Custom object to be saved out.
    dat_type : {'meg', 'maps'}
        Data type of the object.
    save_name : str
        Name to attach to saved out file name.
    db : OMDB() object, optional
        Database object for the OM project.
    """

    # Check db, initialize if not provided
    db = check_db(db)

    # Check that specified dat type is vale
    if dat_type not in ['meg', 'maps']:
        raise UnknownDataTypeError('Data type not understood.')

    # Set file name to save out
    save_name = dat_type + '_' + save_name + '_' + datetime.datetime.now().strftime("%Y-%m-%d") + '.p'

    # Save out data to pickle file
    with open(os.path.join(db.save_path, dat_type, save_name), 'wb') as pickle_file:
        pickle.dump(obj, pickle_file)


def load_obj_pickle(dat_type, file_name, db=None):
    """Load a custom object from a pickled file.

    Parameters
    ----------
    dat_type : {'meg', 'maps'}
        Data type of the object.
    file_name : str
        Label in the filename to be loaded.
    db : OMDB object, optional
        Database object for the OM project.

    Returns
    -------
    pickled object
        Custom object loaded from pickle.
    """

    # Check db, initialize if not provided
    db = check_db(db)

    # Check that specified data type is valid
    if dat_type not in ['meg', 'maps']:
        raise UnknownDataTypeError('Data type not understood.')

    # Check what files are available
    files = os.listdir(os.path.join(db.save_path, dat_type))
    f_names = clean_file_list(files, file_name)

    # Check if there is a single file meeting description
    if len(f_names) == 0:
        raise UnknownDataSourceError('No files found matching description.')
    elif len(f_names) > 1:
        raise UnknownDataSourceError('Multiple files found, be more specific.')
    else:
        f_name = f_names[0]

    with open(os.path.join(db.save_path, dat_type, f_name), 'rb') as pickle_file:
        results = pickle.load(pickle_file)

    return results


def load_meg_list(sub_nums, osc_bands_vert=False, all_oscs=False, osc=None, db=None, dat_source=None):
    """Loads a group of subject IDs into MegSubj() objects, collected in a list.

    Parameters
    ----------
    sub_nums : list of int
        Subject IDs to load and add to list.
    osc_bands_vert : boolean, optional (default: False)
        Whether convert data to oscillations bands across vertices.
    all_oscs : boolean, optional (default: False)
        Whether to convert data to all oscillations.
    osc : Osc() object, optional
        Oscillation band definitions.
    db : OMDB() object, optional
        OM database object.
    dat_source : str, optional
        Which database subject comes from.

    Returns
    -------
    dat_out : list of MegSubj() objects
        A list containing all the loaded MegSubj objects.
    """

    # Check db, initialize if not provided
    db = check_db(db)

    # Loop through subject numbers to load
    dat_out = []
    for subj_id in sub_nums:

        # Initialize MegSubj, load fooof dat
        temp = MegSubj(db, dat_source, osc)
        temp.import_fooof(subj_id, get_demo=True)

        # Convert to oscillation vertex bands, if requested
        if osc_bands_vert:
            temp.osc_bands_vertex()

        # Convert to all-oscs, if requested
        if all_oscs:
            temp.all_oscs(verbose=False)

        # Add new subject to output list
        dat_out.append(temp)

    return dat_out


def load_meg_pairs(dat_list, osc_bands_vert=False, all_oscs=False, osc=None, db=None):
    """Loads pairs of subject IDs into MegSubj() objects.

    Parameters
    ----------
    pairs_list : list of tuple of (int, int)
        Pairs of subject IDs to load and add to list.
    osc_bands_vert : boolean, optional (default: False)
        Whether convert data to oscillations bands across vertices.
    all_oscs : boolean, optional (default: False)
        Whether to convert data to all oscillations.
    osc : Osc() object, optional
        Oscillation band definitions.
    db : OMDB() object, optional
        OM database object.

    Returns
    -------
    list of list of (MegSubj, MegSubj)
        Data for all subject pairs.
    """

    return [load_meg_list(pair, osc_bands_vert=osc_bands_vert, all_oscs=all_oscs,
                          osc=osc, db=db, dat_source='HCP') for pair in dat_list]
