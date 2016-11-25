"""MODULE DOC STRING

"""

from __future__ import print_function, division

import os

import om.cl.md_sing as md_sing
import om.cl.md_gr as md_gr
from om.gen import Osc

#############################################################################################
#############################################################################################
#############################################################################################

class TestDB(object):
    """   """

    def __init__(self):
        """    """

        #
        self.dat_source = 'test'

        #
        base_path = ("/Users/thomasdonoghue/Documents/GitCode/omegamappin/om/tests/data/")
        self.external_path = os.path.join(base_path, 'External')
        self.internal_path = os.path.join(base_path, 'Internal')
        self.other_path = os.path.join(base_path, 'Other')

        # Set up external test data links
        self.meg_path = os.path.join(self.external_path, 'MEG')
        self.foof_path = os.path.join(self.meg_path, 'FOOF')
        self.psd_path = os.path.join(self.meg_path, 'PSDs')
        self.viz_path = os.path.join(self.meg_path, 'Viz')

        # Set up internal test data links
        self.maps_path = os.path.join(self.internal_path, 'Maps')
        self.corrs_path = os.path.join(self.internal_path, 'Corrs')

        # Set
        self.csvs_path = os.path.join(self.other_path, 'csvs')

#############################################################################################
#############################################################################################
#############################################################################################

def load_test_meg_subj(sub):
    """Loads a test subject of MD_SING data."""

    tdb = TestDB()
    osc = Osc(default=True)

    dat = md_sing.MegData(tdb, '', osc)

    dat.import_foof(sub, get_demo=False, load_type='pickle')

    return dat

def load_test_meg_gr(bands_vertex=False, all_osc=False, peaks=False, calc_maps=False, vertex_osc=False):
    """Loads a test group object of MD_GR data."""

    tdb = TestDB()
    osc = Osc(default=True)

    meg_group = md_gr.GroupMegData(tdb, osc)

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
