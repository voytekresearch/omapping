"""   """

from __future__ import print_function, division

import os
import time
from ipyparallel import Client

import sys

# TODO: Add check for 'cluster already started' in file, otherwise can hang there'

####################################################################################
############################# OMEGAMAPPIN - CORE - PAR #############################
####################################################################################

class Par(object):
    """

    Attributes
    ----------
    active : boolean
        xx
    f_name : str
        xx
    verbose : True
        xx
    client : xx
        xx
    workers : xx
        xx
    """

    def __init__(self):
        """   """

        self.active = False
        self.f_name = 'cluster.txt'
        self.verbose = True

        self.client = None
        self.workers = None


    def launch(self, n_core=4):
        """

        Parameters
        ----------
        n_core : int
            Number of cores to run in parallel.
        """

        command = "ipcluster start --n=" + str(n_core) + " &> " + self.f_name + " &"

        os.system(command)

        self.active = True
        time.sleep(0.25)

        self.wait(5, 0.25)

        self.client = Client()
        self.workers = self.client[:]

        if self.verbose:
            print('Cluster opened')


    def stop(self):
        """   """

        os.system("ipcluster stop")

        self.wait(8, 0.1)

        self.active = False

        os.remove(self.f_name)

        if self.verbose:
            print('Cluster shut down.')


    def wait(self, n_lines, n_wait):
        """

        Parameters
        ----------
        n_lines : int
            Number of lines to wait to be in the parallel log file.
        n_wait : float
            Number of seconds to wait each iteration if status is not ready yet.
        """

        while True:

            with open(self.f_name) as f:
                num_lines = sum(1 for line in f)

            if num_lines == n_lines:
                break

            else:
                time.sleep(n_wait)
