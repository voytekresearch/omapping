"""   """

import os
import pickle
import datetime
import numpy as np
import scipy.io as sio

from om.core.utils import clean_file_list

##
##
##


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
        mat_file = 'PSD_Source_Median_' + str(subj_num)
    file_name = os.path.join(meg_path, dat_source, ('Subject_' + str(subj_num)), mat_file)

    # Load MEG PSD data from matfile
    data_mat = sio.loadmat(file_name, appendmat=True, struct_as_record=False, squeeze_me=True)

    # Pull out data from dictionary
    freqs = data_mat['Freqs']
    psd = data_mat['TF']

    return psd, freqs


def save_foof_pickle(results, save_path, sub_num):
    """Save out the FOOF results as a pickle file.

    Parameters
    ----------
    results : list of tuple
        FOOF results - (slope (float), centers (1d array), amps (1d array), bws (1d array)).
    save_path: str
        Filepath of where to save out the file.
    sub_num : int
        Subject identifier number.
    """

    # Set save name and path
    save_name = str(sub_num) + '_Foof_Vertex.p'
    foof_save_path = os.path.join(save_path, 'pickle', save_name)

    # Save out data to pickle file
    pickle.dump(results, open(foof_save_path, 'wb'))


def save_foof_csv(results, save_path, sub_num):
    """Save out the FOOF results as a csv file.

    Parameters
    ----------
    results : list of tuple
        FOOF results - (slope (float), centers (1d array), amps (1d array), bws (1d array)).
    save_path : str
        Filepath of where to save out the file.
    sub_num : int
        Subject identifier number.

    TODO: Update to use csv.
    """

    # Set index values for oscillation parameters
    i_cen = 1
    i_amp = 2
    i_bw = 3

    # Check how many vertices there are
    n_verts = len(results)

    # Initialize file names for oscillation and slope csv files
    csv_sl_fname = save_path + '/csv/' + str(sub_num) + '_Slopes.csv'
    csv_osc_fname = save_path + '/csv/' + str(sub_num) + '_Oscs.csv'

    # Open files to write to
    sl_csv = open(csv_sl_fname, 'w')
    osc_csv = open(csv_osc_fname, 'w')

    # Loop through each vertex
    for vert in range(n_verts):

        # Save out slope value to csv file
        sl_csv.write(str(results[vert][0]) + '\n')

        # Check how oscillations at current vertex
        n_oscs = len(results[vert][1])

        # Loop through each oscillation, writing each one to file
        for osc in range(n_oscs):

            cur_osc_dat = list([vert + 1, results[vert][i_cen][osc],
                                results[vert][i_amp][osc], results[vert][i_bw][osc]])
            osc_csv.write((", ".join(repr(el) for el in cur_osc_dat)) + '\n')


def load_foof_pickle(dat_path, sub_num):
    """Load FOOF data from pickle file.

    Parameters
    ----------
    dat_path : str
        File name for where data is stored to load from.
    sub_num : int
        Subject identifier number.

    Returns
    -------
    results : list of tuple
        FOOF results - (slope (float), centers (1d array), amps (1d array), bws (1d array)).
    """

    # Get list of available files to load
    files = os.listdir(os.path.join(dat_path, 'pickle'))
    files = clean_file_list(files, 'Foof_Vertex')

    # Get specific subject file
    cur_subj_file = get_cur_subj(sub_num, files)
    subj_path = os.path.join(dat_path, 'pickle', cur_subj_file)

    # Load file
    results = pickle.load(open(subj_path, 'rb'))

    return results


def load_foof_csv():
    """
    NOTE: Not yet implemented.
    """

    pass

def save_md_pickle(obj, save_name, db=None):
    """Save current meg data object as a pickled object.

    Parameters
    ----------
    obj : MegData() or GroupMegData()
        Object to save to pickle
    save_name : str
        String to be included in the name of the file.
    """

    # Get database object, unless one was provided
    if not db:
        db = OMDB()

    # Set save name and path
    save_name = 'Res_' + save_name + '_' + datetime.datetime.now().strftime("%Y-%m-%d") + '.p'

    # Save out data to pickle file
    pickle.dump(obj, open(os.path.join(db.md_save_path, save_name), 'wb'))

def load_md_pickle(file_name, db=None):
    """Load a pickled file.

    Parameters
    ----------
    file_name : str
        File name of the pickle file to be loaded.

    Returns
    -------
    results : ?
        xx
    """

    # Get database object, unless one was provided
    if not db:
        db = OMDB()

    # Check what files are available
    files = os.listdir(db.md_save_path)
    f_names = clean_file_list(files, file_name)

    # Check if there is a single file meeting description
    if len(f_names) > 1:
        raise UnknownDataSourceError('Unclear which file to load - be more specific.')
    else:
        f_name = f_names[0]

    # Load file & return pickled object
    return pickle.load(open(os.path.join(db.md_save_path, f_name), 'rb'))
