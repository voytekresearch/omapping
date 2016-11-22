"""   """

from __future__ import print_function, division

import os

import om.cl.md_sing as md
from om.gen import Osc

##
##
##

class TestDB(object):

    def __init__(self):
        """    """

        self.dat_source = 'test'
        self.project_path = ("/Users/thomasdonoghue/Documents/GitCode/omegamappin/om/tests/test_files/")

        self.foof_path = os.path.join(self.project_path, 'foof')

##
##
##

def load_test_meg_subj(sub):
    """Loads a test subject of MD_SING data."""

    db = TestDB()
    osc = Osc(default=True)

    dat = md.MegData(db, osc)

    dat.import_foof(sub, get_demo=False, load_type='pickle')

    return dat