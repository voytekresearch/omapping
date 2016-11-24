"""DOCSTRING - TO FILL IN

"""

from __future__ import print_function, division
import os
import sys
import time
import numpy as np

# Import required things from ipyparallel
from ipyparallel import Client
from ipyparallel.util import interactive

# Import general code from custom module om
sys.path.append('/Users/thomasdonoghue/Documents/GitCode/omegamappin/')
from om.gen import OMDB
from om.gen import clean_file_list, get_sub_nums, load_meg_psds, save_foof_pickle, extract_psd

## TODO:
# - Add a report to save out
# - Figure out ipyparallel to launch from within script
# - Turn into proper python script with main() function

############################################################################
################################# SETTINGS #################################
############################################################################

# Set which data set to run. Options {'HCP', 'OMEGA'}
DAT_SOURCE = 'HCP'

# Initiate queue of subjects to process
MEG_QUEUE = [358144, 433839, 512835, 555348, 559053, 568963, 581450, 599671]

###########################################################################
################################ FUNCTIONS ################################
###########################################################################

# Define function to run foof
@interactive
def run_foof(psd_in):
    """DOCSTRING"""

    # Initialize foof object
    foof = FOOF(min_p=min_p, res=freq_res, fmin=fmin, fmax=fmax)

    # Model foof
    foof.model(freqs_ext, psd_in)

    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)

############################################################################
################################# RUN CODE #################################
############################################################################

def main():

    # Print out beginning status
    print('\n\nSTARTING FOOF ON MEG DATA: \n')
    time.sleep(2)

    ## Set up parallel workers
    # Initialize client
    c = Client()
    # Gather
    view = c[:]
    print('Workers Connected... \n')

    # Import required libraries to each worker
    print('Doing engine imports...')
    with view.sync_imports():
        import sys
        sys.path.append('/Users/thomasdonoghue/Documents/GitCode/')
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
        psd, freqs = load_meg_psds(db.DAT_SOURCE, db.psd_path, subj)

        # Check data - get number of PSDs and frequency resolution
        [nPSDs, nFreqs] = np.shape(psd)
        freq_res = np.mean(np.diff(freqs))

        # Extract PSD range of interest
        psd_ext, freqs_ext = extract_psd(psd, freqs, f_low=3, f_high=40)

        # Set up PSD as a list of 2-D np arrays
        psd_list = list(psd_ext)
        for i in range(nPSDs):
            psd_list[i] = np.reshape(psd_list[i], [len(freqs_ext), 1])

        # Send required vars to workers
        view['min_p'] = 0.1
        view['freq_res'] = freq_res
        view['fmin'] = freqs_ext.min()
        view['fmax'] = freqs_ext.max()
        view['freqs_ext'] = freqs_ext

        # Set up and run foof parallel
        foof_map = view.map(run_foof, psd_list)
        foof_results = foof_map.get()

        # Save out results
        save_foof_pickle(foof_results, os.path.join(db.foof_path, DAT_SOURCE), subj)

        # Print status
        print('FOOF finished and saved for subj ', str(subj))
        print('Finished at ', time.strftime('%H:%M:%S', time.localtime()), '\n')

        # Pause before next subject
        time.sleep(200)

    # Print out end status
    print('\nFINISHED FOOF ON MEG DATA\n\n')


if __name__ == "__main__":
    main()
