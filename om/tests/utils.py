"""Utility and helper functions for the testing of OM."""

import os

#import om.cl.md_sing as md_sing
from om.meg.single import MegData
#import om.cl.md_gr as md_gr
from om.meg.group import GroupMegData

#from om.gen import OMDB, Osc
from om.core.db import OMDB
from om.core.osc import Osc

#############################################################################################
#############################################################################################
#############################################################################################

class TestDB(OMDB):
    """   """

    def __init__(self):

        # Initialize from OMDB object
        OMDB.__init__(self)

        # Reset base paths
        self.internal_path = ("/Users/thomasdonoghue/Documents/GitCode/"
                              "omegamappin/om/tests/data/Internal/")
        self.external_path = ("/Users/thomasdonoghue/Documents/GitCode/"
                              "omegamappin/om/tests/data/External/")

        # Generate test paths
        self.gen_paths()

        # Add other paths specific to testing
        self.other_path = ("/Users/thomasdonoghue/Documents/GitCode/"
                           "omegamappin/om/tests/data/Other/")
        self.csvs_path = os.path.join(self.other_path, 'csvs')

#############################################################################################
#############################################################################################
#############################################################################################

def load_test_meg_subj(sub):
    """Loads a test subject of MD_SING data."""

    tdb = TestDB()
    osc = Osc(default=True)

    dat = MegData(tdb, '', osc)

    dat.import_foof(sub, get_demo=False, load_type='pickle')

    return dat

def load_test_meg_gr(bands_vertex=False, all_osc=False, peaks=False, calc_maps=False, vertex_osc=False):
    """Loads a test group object of MD_GR data."""

    tdb = TestDB()
    osc = Osc(default=True)

    meg_group = GroupMegData(tdb, osc)

    subjs = ['test_v5', 'test_v5']

    for s in subjs:

        meg_subj = load_test_meg_subj(s)

        if bands_vertex:
            meg_subj.osc_bands_vertex()

        if all_osc:
            meg_subj.all_oscs()

        if peaks:
            meg_subj.peak_freq()

        meg_group.add_subject(meg_subj, add_all_oscs=all_osc,
                              add_vertex_bands=bands_vertex,
                              add_peak_freqs=peaks,
                              add_vertex_oscs=vertex_osc)

    if calc_maps:
        meg_group.osc_prob()
        meg_group.osc_score()

    return meg_group
