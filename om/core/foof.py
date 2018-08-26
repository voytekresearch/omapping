"""Functions to run FOOF for the OM project.

NOTE: This code is old (unused - totally outdated.)
"""

import sys
import numpy as np

# Import FOOF (use sys to add location to path, then import)
#sys.path.append('/Users/tom/Documents/GitCode/')
#from foof.fit import FOOF

###################################################################################
################################ OM_FOOF Functions ################################
###################################################################################

def meg_foof(psd_ext, freqs_ext, min_p, freq_res):
    """Run FOOF on MEG-PSD data.

    NOTE: This is an old version, can run foof linearly.
        For parallel, use stand-alone script.

    Parameters
    ----------
    psd_ext : 2d array
        Matrix of PSDs in the form of [n_verts, n_freqs].
    freqs_ext : 1d array
        Vector of the frequency values for each power value in psd_ext.
    min_p : float
        Minimum probability for splitting peaks. Parameter for FOOF.
    freqs_res : float
        Frequency resolution.

    Returns
    -------
    results : list of tuple
        FOOF results - (slope (float), centers (1d array), amps (1d array), bws (1d array)).
    """

    # Check how many PSDs there are
    [n_PSDs, n_freqs] = np.shape(psd_ext)

    # Initialize foof
    foof = FOOF(min_p=min_p, res=freq_res, fmin=freqs_ext.min(), fmax=freqs_ext.max())

    # Set up PSD as a list of 2-D np arrays
    psd_list = list(psd_ext)
    for i in range(n_PSDs):
        psd_list[i] = np.reshape(psd_list[i], [len(freqs_ext), 1])

    # Run FOOF linearly
    results = [_run_foof_l(foof, freqs_ext, psd) for psd in psd_list]

    return results

########################################################################################
############################## OM GEN - Private Functions ##############################
########################################################################################

def _run_foof_l(foof, freqs_ext, psd_ext):
    """Local helper function to run FOOF linearly.

    Used by meg_foof().

    Parameters
    ----------
    foof : FOOF() object
        FOOF object to model 1/f & oscillations.
    freqs_ext : 1d array
        Vector of frequency values for the psd.
    psd_ext : 1d array
        Vector of power values for the psd.

    Returns
    -------
    out : tuple
        FOOF results - (slope (float), centers (1d array), amps (1d array), bws (1d array)).
    """

    # Fit FOOF
    foof.model(freqs_ext, psd_ext)

    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)
