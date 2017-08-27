"""Parallel functions for OM project."""

from __future__ import print_function, division

import os
import time
from ipyparallel import Client
from ipyparallel.util import interactive

from om.core.errors import ClusterAlreadyRunningError

####################################################################################
############################# OMEGAMAPPIN - CORE - PAR #############################
####################################################################################

class Par(object):
    """Object to control and support parallel computations.

    Attributes
    ----------
    active : boolean
        Whether or not the parallel cluster is active.
    f_name : str
        Name of the log file.
    verbose : True
        Whether to print out status reports.
    client : Client()
        Parallel client.
    workers : ?
        Parallel worker engines.
    """

    def __init__(self, verbose=True):
        """Intialize parallel object."""

        self.active = False
        self.f_name = 'cluster.txt'
        self.verbose = verbose

        self.client = None
        self.workers = None


    def launch(self, n_core=4):
        """Initiate parallel workers.

        Parameters
        ----------
        n_core : int
            Number of cores to run in parallel.
        """

        if self.verbose:
            print('\n Starting Cluster...')

        # Send command to start cluster
        command = "ipcluster start --n=" + str(n_core) + " &> " + self.f_name + " &"
        os.system(command)
        time.sleep(0.25)

        # Check if a cluster is already open
        self.check_for_open()

        # Set as active, wait for it to finish
        self.active = True
        self.wait(5, 0.25)

        # Collect worker engines
        self.client = Client()
        self.workers = self.client[:]

        if self.verbose:
            print('Cluster Started. \n')


    def stop(self):
        """Stop parallel workers."""

        if self.verbose:
            print('\n Shutting down Cluster...')

        # Send command to stop cluster, wait for it to finish
        os.system("ipcluster stop")
        self.wait(8, 0.1)

        # Remove the log file, and set status
        os.remove(self.f_name)
        self.active = False

        if self.verbose:
            print('Cluster shut down. \n')


    def wait(self, n_lines, n_wait):
        """Wait for log file.

        Parameters
        ----------
        n_lines : int
            Number of lines to wait to be in the parallel log file.
        n_wait : float
            Number of seconds to wait each iteration if status is not ready yet.
        """

        while True:

            # Check status of log file
            with open(self.f_name) as cur_file:
                num_lines = sum(1 for line in cur_file)

            # Check if log has required number of lines
            if num_lines == n_lines:
                break

            # Otherwise, wait and re-try
            else:
                self.check_for_open()
                time.sleep(n_wait)


    def check_for_open(self):
        """Check if a cluster is already open, while trying to open a new one."""

        # Use log file to check for text containing info that a cluster is already up
        with open(self.f_name) as cur_file:

            # If any row of the log file has the message, throw an error
            for row in cur_file:
                if 'Cluster is already running' in row:
                    raise ClusterAlreadyRunningError("Can't start cluster, already running.")

###################################################################################
########################## OMEGAMAPPIN - PAR - FUNCTIONS ##########################
###################################################################################

@interactive
def run_foof_par(psd_in):
    """Function to pass to run FOOF in parallel.

    Parameters
    ----------
    psd_in : ?
        xx

    Notes
    -----
    - FOOF has to be imported on workers.
    - min_p, freq_res, fmin, fmax have to be sent to workers.

    NOTE: FOR OLD FOOF
    """

    # Initialize foof object
    foof = FOOF(min_p=min_p, res=freq_res, fmin=fmin, fmax=fmax)

    # Model foof
    foof.model(freqs_ext, psd_in)

    # Store vals in tuple and return
    return (foof.chi_, foof.centers_, foof.powers_, foof.stdevs_)


@interactive
def run_fooof_par(psd_in):
    """Function to pass to run FOOOF in parallel.

    Parameters
    ----------
    psd_in : 1d array
        PSD to run FOOOF on.

    Returns
    -------
    tuple
        FOOOF results.

    Notes
    -----
    - FOOOF has to be imported on workers.
    - freqs, freq_range have to be sent to workers.
    """

    # Initialize FOOOF object
    fm = FOOOF(bandwidth_limits=bandwidth_limits, max_n_oscs=max_n_oscs)

    # Fit the PSD model
    fm.fit(freqs, psd_in, freq_range=freq_range)

    # Return FOOOF model fit parameters
    return (fm.get_params())


@interactive
def run_corr_par(dat):
    """Run correlation between maps. Used for parallel runs.

    Parameters
    ----------
    dat : 1d array
        An array of map data to be compared to projected meg map.

    Returns
    -------
    out : tuple
        Correlation results (corr_vals, p_vals).

    Notes:
    - meg_map has to be projected to workers.
    - numpy and pearsonr have to be imported on workers.
    """

    # Get inds of data that contains numbers
    inds_non_nan = numpy.invert(numpy.isnan(dat))

    # Calculate corr between data and MEG map
    [corr_vals, p_vals] = pearsonr(dat[inds_non_nan], meg_map[inds_non_nan])

    return (corr_vals, p_vals)
