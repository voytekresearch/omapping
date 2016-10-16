from __future__ import print_function, division

# Import required packages
import os
import sys
import time
import numpy as np

# Import general code from custom module om
sys.path.append('/Users/thomasdonoghue/Documents/GitCode/omegamappin/')
from om.gen import *

# Import required things from ipyparallel
from ipyparallel import Client
from ipyparallel.util import interactive

## TODO: Add a report to save out.

######################################################################
############################## SETTINGS ##############################
######################################################################

# Set which data set to run. Options {'HCP', 'OMEGA'}
dat_source = 'OMEGA'

# Initiate queue of subjects to process
meg_queue = [559176, 369737]

###################################################################
############################ FUNCTIONS ############################
###################################################################

# Define function to run foof
@interactive
def run_foof(psd_in):

    # Initialize foof object
    foof = FOOF(min_p=min_p, res=freq_res, fmin=fmin, fmax=fmax)

    # Model foof
    foof.model(freqs_ext, psd_in)

    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)

############################################################################
################################# RUN CODE #################################
############################################################################

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

## Set Paths to MEG Data
# OMEGA
if dat_source is 'OMEGA':
    meg_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/2-PSDs/OMEGA/'
    save_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/3-FOOF/OMEGA/pickle'
# HCP
elif dat_source is 'HCP':
    meg_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/2-PSDs/HCP/'
    save_path = '/Users/thomasdonoghue/Documents/Research/1-Projects/OMEGA/2-Data/MEG/3-FOOF/HCP/pickle'
else:
    print('Data source not understood.')
    time.sleep(5)

# Check Availabe Subjects
files = os.listdir(meg_path)
files = clean_file_list(files, 'Subject_')
all_sub_nums = get_sub_nums(files, 'last')

# Get list of subjects who have already been FOOFed
foofed_subj_files = os.listdir(save_path)
foofed_subj_files = clean_file_list(foofed_subj_files, 'foof')
foofed_subj_nums = get_sub_nums(foofed_subj_files, 'first')

# Check all subject numbers given are unique
if len(set(all_sub_nums)) != len(all_sub_nums):
    print('There is a duplicate in the given subj numbers!')
    time.sleep(2)

# Check each subject is avaialble, and whether it has already been processed
for sub_num in meg_queue:

    # Check if subject is available
    if sub_num not in all_sub_nums:
        print('\nSubject ', str(sub_num), ' not found. Removing from list.')
        meg_queue.remove(sub_num)

    # Check if subject has already been FOOFed
    if sub_num in foofed_subj_nums:
        print('\nSubject ', str(sub_num), ' has already been FOOFed.')

# Print out subject status
print('\nSet to run subjects: ', meg_queue, '\n')

# Loop through subjects set to run
for subj in meg_queue:

    # Print out status
    print('Running FOOF on subj ', str(subj))
    print('Starting at ', time.strftime('%H:%M:%S', time.localtime()))

    # Load MEG Data
    psd, freqs = load_meg_psds(meg_path, subj)

    # Check data - get number of PSDs and frequency resolution
    [nPSDs, nFreqs] = np.shape(psd)
    freq_res = np.mean(np.diff(freqs))

    # Extract PSD range of interest
    psd_ext, freqs_ext = extract_psd(psd, freqs, f_low=3, f_high=40)

    # Set up PSD as a list of 2-D np arrays
    psd_list = list(psd_ext)
    for i in range(0, nPSDs):
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
    save_foof_pickle(foof_results, save_path, subj)

    # Print status
    print('FOOF finished and saved for subj ', str(subj))
    print('Finished at ', time.strftime('%H:%M:%S', time.localtime()), '\n')

    # Pause before next subject
    time.sleep(200)

# Print out end status
print('\nFINISHED FOOF ON MEG DATA\n\n')


