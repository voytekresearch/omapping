"""   """

from __future__ import print_function, division

import os

import om.cl.md_sing as md_sing
import om.cl.md_gr as md_gr
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

        self.csvs_path = os.path.join(self.project_path, 'csvs')

        self.maps_path = os.path.join(self.project_path, 'Maps')

##
##
##

def load_test_meg_subj(sub):
    """Loads a test subject of MD_SING data."""

    tdb = TestDB()
    osc = Osc(default=True)

    dat = md_sing.MegData(tdb, '', osc)

    dat.import_foof(sub, get_demo=False, load_type='pickle')

    return dat

def load_test_meg_gr(bands_vertex=False):
    """Loads a test group object of MD_GR data."""

    tdb = TestDB()
    osc = Osc(default=True)

    meg_group = md_gr.GroupMegData(tdb, osc)

    subjs = ['test2', 'test2']

    for s in subjs:

        meg_subj = load_test_meg_subj(s)

        if bands_vertex:
            meg_subj.osc_bands_vertex()

        meg_group.add_subject(meg_subj, add_vertex_bands=bands_vertex)

    return meg_group
