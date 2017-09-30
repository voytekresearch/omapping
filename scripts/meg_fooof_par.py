"""Standalone script for running FOOOF on a group of MEG subjects (in parallel)."""

import os
import time
import numpy as np

from om.core.par import Par, run_fooof_par
from om.core.db import OMDB
from om.core.io import load_meg_psds, save_fooof_pickle
from om.core.utils import clean_file_list, get_sub_nums

######################################################################################
###################################### SETTINGS ######################################
######################################################################################

# Set which data set to run. Options {'HCP', 'OMEGA'}
DAT_SOURCE = 'HCP'

# Initiate queue of subjects to process
MEG_QUEUE = [191437, 191841]
            # , 192641, 195041, 198653,
            #  204521, 205119, 212318, 212823, 214524, 221319, 223929, 233326, 248339,
            #  250427, 255639, 257845, 283543, 287248, 293748, 352132, 352738, 353740,
            #  358144, 406836, 433839, 512835, 555348, 559053, 568963, 581450, 599671,
            #  601127, 660951, 662551, 665254, 667056, 679770, 706040, 707749, 715950,
            #  725751, 735148, 783462, 814649, 825048, 877168, 891667, 898176, 912447,
            #  917255, 990366]

# FOOOF settings
FREQ_RANGE = [3, 40]
BANDWIDTH_LIMITS = [1, 8]
MAX_N_OSCS = 8
FIT_KNEE = True

######################################################################################
###################################### RUN CODE ######################################
######################################################################################

def main():
    """Run FOOOF in parallel across MEG subjects."""

    # Print out beginning status
    print('\n\nSTARTING FOOOF ON MEG DATA: \n')
    time.sleep(2)

    ## Set up parallel workers
    par = Par()
    par.launch()

    # Import required libraries to each worker
    print('Doing engine imports...')
    with par.workers.sync_imports():
        from fooof import FOOOF

    # Set up database object
    db = OMDB()

    # Check Availabe Subjects
    files = os.listdir(os.path.join(db.psd_path, DAT_SOURCE))
    files = clean_file_list(files, 'Subject_')
    all_sub_nums = get_sub_nums(files, 'last')

    # Get list of subjects who have already been FOOOFed
    fooofed_subj_files = os.listdir(os.path.join(db.fooof_path, DAT_SOURCE, 'pickle'))
    fooofed_subj_files = clean_file_list(fooofed_subj_files, 'fooof')
    fooofed_subj_nums = get_sub_nums(fooofed_subj_files, 'first')

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

        # Check if subject has already been fooofed
        if sub_num in fooofed_subj_nums:
            print('\nSubject ', str(sub_num), ' has already been FOOOFed.')

    # Print out subject status
    print('\nSet to run subjects: ', MEG_QUEUE, '\n')

    # Loop through subjects set to run
    for subj in MEG_QUEUE:

        # Print out status
        print('Running FOOOF on subj ', str(subj))
        print('Starting at ', time.strftime('%H:%M:%S', time.localtime()))

        # Load MEG Data
        psd, freqs = load_meg_psds(DAT_SOURCE, db.psd_path, subj)

        # Set up PSD as a list of 1d np arrays
        psd_list = [cur_psd for cur_psd in psd]

        # Send required vars to workers
        par.workers['freqs'] = freqs
        par.workers['freq_range'] = FREQ_RANGE
        par.workers['bandwidth_limits'] = BANDWIDTH_LIMITS
        par.workers['max_n_oscs'] = MAX_N_OSCS
        par.workers['fit_knee'] = FIT_KNEE

        # Set up and run fooof parallel
        fooof_map = par.workers.map(run_fooof_par, psd_list)
        fooof_results = fooof_map.get()

        # Save out results
        save_fooof_pickle(fooof_results, os.path.join(db.fooof_path, DAT_SOURCE), subj)

        # Print status
        print('fooof finished and saved for subj ', str(subj))
        print('Finished at ', time.strftime('%H:%M:%S', time.localtime()), '\n')

        # Pause before next subject
        time.sleep(180)

    # Print out end status
    par.stop()
    print('\nFINISHED fooof ON MEG DATA\n\n')


if __name__ == "__main__":
    main()
