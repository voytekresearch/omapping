"""DOCSTRING - TO FILL IN"""

from __future__ import print_function, division

import os
import sys
import time
import numpy as np

# Import required things from ipyparallel
from ipyparallel import Client
from ipyparallel.util import interactive

# Import general code from custom module om
from om.core.par import Par, run_foof_par
from om.core.db import OMDB
from om.core.io import load_meg_psds, save_foof_pickle
from om.core.utils import clean_file_list, get_sub_nums, extract_psd

############################################################################
################################# SETTINGS #################################
############################################################################

# Set which data set to run. Options {'HCP', 'OMEGA'}
DAT_SOURCE = 'HCP'

# Initiate queue of subjects to process
MEG_QUEUE = [358144, 433839, 512835, 555348, 559053, 568963, 581450, 599671]

############################################################################
################################# RUN CODE #################################
############################################################################

def main():
    """DOCSTRING"""

    # Print out beginning status
    print('\n\nSTARTING FOOF ON MEG DATA: \n')
    time.sleep(2)

    ## Set up parallel workers
    par = Par()
    par.launch()

    # Import required libraries to each worker
    print('Doing engine imports...')
    with par.workers.sync_imports():
        from foof.fit import FOOF

    # Set up database object
    db = OMDB()

    # Check Availabe Subjects
    files = os.listdir(os.path.join(db.psd_path, DAT_SOURCE))
    files = clean_file_list(files, 'Subject_')
    all_sub_nums = get_sub_nums(files, 'last')

    # Get list of subjects who have already been FOOFed
    foofed_subj_files = os.listdir(os.path.join(db.foof_path, DAT_SOURCE, 'pickle'))
    foofed_subj_files = clean_file_list(foofed_subj_files, 'foof')
    foofed_subj_nums = get_sub_nums(foofed_subj_files, 'first')

    # Check all subject numbers given are unique
    if len(set(all_sub_nums)) != len(all_sub_nums):
        print('There is a duplicate in the given subj numbers!')
        time.sleep(2)

    # Check each subject is avaialble, and whether it has already been processed
    for sub_num in MEG_QUEUE:

        # Check if subject is available
        if sub_num not in all_sub_nums:
            print('\nSubject ', str(sub_num), ' not found. Removing from list.')
            MEG_QUEUE.remove(sub_num)

        # Check if subject has already been FOOFed
        if sub_num in foofed_subj_nums:
            print('\nSubject ', str(sub_num), ' has already been FOOFed.')

    # Print out subject status
    print('\nSet to run subjects: ', MEG_QUEUE, '\n')

    # Loop through subjects set to run
    for subj in MEG_QUEUE:

        # Print out status
        print('Running FOOF on subj ', str(subj))
        print('Starting at ', time.strftime('%H:%M:%S', time.localtime()))

        # Load MEG Data
        psd, freqs = load_meg_psds(DAT_SOURCE, db.psd_path, subj)

        # Check data - get number of PSDs and frequency resolution
        [n_psds, n_freqs] = np.shape(psd)
        freq_res = np.mean(np.diff(freqs))

        # Extract PSD range of interest
        psd_ext, freqs_ext = extract_psd(psd, freqs, f_low=3, f_high=40)

        # Set up PSD as a list of 2-D np arrays
        psd_list = list(psd_ext)
        for i in range(n_psds):
            psd_list[i] = np.reshape(psd_list[i], [len(freqs_ext), 1])

        # Send required vars to workers
        par.workers['min_p'] = 0.1
        par.workers['freq_res'] = freq_res
        par.workers['fmin'] = freqs_ext.min()
        par.workers['fmax'] = freqs_ext.max()
        par.workers['freqs_ext'] = freqs_ext

        # Set up and run foof parallel
        foof_map = par.workers.map(run_foof_par, psd_list)
        foof_results = foof_map.get()

        # Save out results
        save_foof_pickle(foof_results, os.path.join(db.foof_path, DAT_SOURCE), subj)

        # Print status
        print('FOOF finished and saved for subj ', str(subj))
        print('Finished at ', time.strftime('%H:%M:%S', time.localtime()), '\n')

        # Pause before next subject
        time.sleep(200)

    # Print out end status
    par.stop()
    print('\nFINISHED FOOF ON MEG DATA\n\n')


if __name__ == "__main__":
    main()
