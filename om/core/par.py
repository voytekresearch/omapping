"""DOCSTRING"""

from __future__ import print_function, division

import os
import time
from ipyparallel import Client

# TODO: Add check for 'cluster already started' in file, otherwise can hang there'

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
    workers : xx
        Parallel worker engines.
    """

    def __init__(self):
        """Intialize parallel object."""

        self.active = False
        self.f_name = 'cluster.txt'
        self.verbose = True

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
                time.sleep(n_wait)
